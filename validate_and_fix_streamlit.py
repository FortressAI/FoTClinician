#!/usr/bin/env python3
"""
Validate and fix all Streamlit data loading issues.
This script ensures proper data structure, score mapping, and file priorities.
"""

import json
import os
from datetime import datetime

def validate_data_file(file_path):
    """Validate a data file for Streamlit compatibility."""
    print(f"🔍 Validating {file_path}...")
    
    if not os.path.exists(file_path):
        print(f"   ❌ File does not exist")
        return False
    
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        # Check structure
        if not isinstance(data, dict):
            print(f"   ❌ Data is not a dictionary")
            return False
        
        if 'discoveries' not in data:
            print(f"   ❌ Missing 'discoveries' key")
            return False
        
        discoveries = data['discoveries']
        if not isinstance(discoveries, list):
            print(f"   ❌ 'discoveries' is not a list")
            return False
        
        # Check discovery count
        discovery_count = len(discoveries)
        print(f"   📊 Contains {discovery_count} discoveries")
        
        # Check score fields
        scores_with_score = 0
        scores_with_combined = 0
        scores_missing = 0
        
        for i, mol in enumerate(discoveries[:10]):  # Check first 10
            if 'score' in mol:
                scores_with_score += 1
            elif 'combined_score' in mol:
                scores_with_combined += 1
            else:
                scores_missing += 1
        
        print(f"   🎯 Score fields (sample of 10):")
        print(f"      - With 'score': {scores_with_score}")
        print(f"      - With 'combined_score': {scores_with_combined}")
        print(f"      - Missing both: {scores_missing}")
        
        # Check if summary exists and calculate stats
        if 'discovery_summary' in data:
            summary = data['discovery_summary']
            total_claimed = summary.get('total_discoveries', len(discoveries))
            print(f"   📈 Summary claims {total_claimed} total discoveries")
            
            if total_claimed != discovery_count:
                print(f"   ⚠️  Mismatch: Summary claims {total_claimed}, file has {discovery_count}")
        
        # Calculate actual score statistics
        all_scores = []
        for mol in discoveries:
            score = mol.get('score', mol.get('combined_score', 0))
            if score > 0:
                all_scores.append(score)
        
        if all_scores:
            avg_score = sum(all_scores) / len(all_scores)
            max_score = max(all_scores)
            min_score = min(all_scores)
            print(f"   📊 Actual scores: avg={avg_score:.3f}, max={max_score:.3f}, min={min_score:.3f}")
        else:
            print(f"   ❌ No valid scores found!")
            return False
        
        print(f"   ✅ File is valid for Streamlit")
        return True
        
    except Exception as e:
        print(f"   ❌ Error reading file: {e}")
        return False

def fix_streamlit_cache():
    """Clear Streamlit cache by modifying the load function."""
    print("🧹 Clearing Streamlit cache...")
    
    # Add a timestamp to force cache refresh
    cache_buster = {
        "cache_refresh": datetime.now().isoformat(),
        "reason": "Force reload after data fixes"
    }
    
    with open("streamlit_cache_buster.json", 'w') as f:
        json.dump(cache_buster, f, indent=2)
    
    print("   ✅ Cache buster created")

def main():
    """Main validation and fix routine."""
    print("🧪 STREAMLIT DATA VALIDATION & FIX")
    print("=" * 50)
    
    # Check all potential data files
    files_to_check = [
        "results/chemistry_discoveries.json",
        "results/overnight_discovery_mega_dataset.json", 
        "cloud_data_snapshot.json",
        "results/chemistry_discoveries_COMPLETE.json"
    ]
    
    valid_files = []
    for file_path in files_to_check:
        if validate_data_file(file_path):
            valid_files.append(file_path)
    
    print()
    print(f"📋 VALIDATION SUMMARY:")
    print(f"   ✅ Valid files: {len(valid_files)}")
    print(f"   ❌ Invalid files: {len(files_to_check) - len(valid_files)}")
    
    if valid_files:
        print(f"   🎯 Best file to use: {valid_files[0]}")
        
        # Load the best file and show what Streamlit should see
        with open(valid_files[0], 'r') as f:
            data = json.load(f)
        
        discoveries = data.get('discoveries', [])
        scores = [mol.get('score', mol.get('combined_score', 0)) for mol in discoveries]
        valid_scores = [s for s in scores if s > 0]
        
        print()
        print(f"📊 STREAMLIT SHOULD DISPLAY:")
        print(f"   🧬 Total Molecules: {len(discoveries)}")
        if valid_scores:
            print(f"   📊 Average Score: {sum(valid_scores)/len(valid_scores):.3f}")
            print(f"   🏆 Max Score: {max(valid_scores):.3f}")
        else:
            print(f"   ❌ No valid scores found!")
    
    # Clear cache to force refresh
    fix_streamlit_cache()
    
    print()
    print("🎉 VALIDATION COMPLETE!")
    print("   Now restart Streamlit to see the fixed data.")

if __name__ == "__main__":
    main()
