#!/usr/bin/env python3
"""
Fix Streamlit display issues:
1. Map combined_score to score for display
2. Export ALL 2,983 discoveries instead of just 100 samples
3. Create proper data structure for Streamlit
"""

import json
import glob
from datetime import datetime
import os

def main():
    print("🔧 FIXING STREAMLIT DISPLAY ISSUES")
    print("=" * 50)
    
    # Load all individual discovery files
    discovery_files = glob.glob("continuous_chemistry_discoveries/discoveries/*.json")
    print(f"📂 Found {len(discovery_files)} individual discovery files")
    
    all_discoveries = []
    total_scores = []
    
    # Process each discovery file
    for i, file_path in enumerate(discovery_files):
        try:
            with open(file_path, 'r') as f:
                discovery = json.load(f)
            
            # Map combined_score to score for Streamlit compatibility
            if 'combined_score' in discovery and 'score' not in discovery:
                discovery['score'] = discovery['combined_score']
            
            # Ensure we have a score
            score = discovery.get('score', discovery.get('combined_score', 0))
            if score:
                total_scores.append(score)
            
            all_discoveries.append(discovery)
            
            if i % 500 == 0:
                print(f"   Processed {i}/{len(discovery_files)} files...")
                
        except Exception as e:
            print(f"⚠️  Error processing {file_path}: {e}")
            continue
    
    print(f"✅ Processed {len(all_discoveries)} discoveries")
    print(f"📊 Score range: {min(total_scores):.3f} - {max(total_scores):.3f}")
    print(f"📊 Average score: {sum(total_scores)/len(total_scores):.3f}")
    
    # Create comprehensive export structure
    export_data = {
        "discovery_summary": {
            "total_discoveries": len(all_discoveries),
            "generated_at": datetime.now().isoformat(),
            "discovery_campaign": "Complete FoTChemistry Discovery Dataset", 
            "statistics": {
                "avg_score": sum(total_scores) / len(total_scores) if total_scores else 0,
                "max_score": max(total_scores) if total_scores else 0,
                "min_score": min(total_scores) if total_scores else 0,
                "unique_molecules": len(set(d.get('smiles', '') for d in all_discoveries)),
                "score_field_mapped": True
            },
            "data_fixes": [
                "Mapped combined_score to score for Streamlit compatibility",
                "Included ALL discoveries (not just sample)",
                "Fixed data structure for proper display"
            ]
        },
        "discoveries": all_discoveries  # ALL discoveries, not just 100
    }
    
    # Save the complete dataset
    output_file = "results/chemistry_discoveries_COMPLETE.json"
    with open(output_file, 'w') as f:
        json.dump(export_data, f, indent=2)
    
    print(f"💾 Saved complete dataset to {output_file}")
    print(f"📦 File size: {os.path.getsize(output_file):,} bytes")
    
    # Also update the main file that Streamlit uses
    main_file = "results/chemistry_discoveries.json"
    with open(main_file, 'w') as f:
        json.dump(export_data, f, indent=2)
    
    print(f"✅ Updated main Streamlit file: {main_file}")
    
    print()
    print("🎉 FIXES COMPLETE!")
    print(f"   📊 Now shows: {len(all_discoveries)} discoveries (was 100)")
    print(f"   🎯 Score display: Fixed (was 0.000)")
    print(f"   📈 Real avg score: {sum(total_scores)/len(total_scores):.3f}")
    
if __name__ == "__main__":
    main()
