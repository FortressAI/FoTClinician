#!/usr/bin/env python3
"""
Update discovery data files with all new discoveries
Consolidates individual discovery files into the main data files
"""

import json
import os
from glob import glob
from datetime import datetime
import sys

def load_individual_discoveries():
    """Load all individual discovery files"""
    discovery_files = glob('continuous_chemistry_discoveries/discoveries/*.json')
    discoveries = []
    
    print(f"📂 Found {len(discovery_files)} individual discovery files")
    
    for i, file_path in enumerate(discovery_files):
        try:
            with open(file_path, 'r') as f:
                discovery = json.load(f)
            discoveries.append(discovery)
            
            if (i + 1) % 1000 == 0:
                print(f"   ✅ Loaded {i + 1}/{len(discovery_files)} discoveries")
                
        except Exception as e:
            print(f"   ⚠️ Failed to load {file_path}: {e}")
    
    print(f"✅ Successfully loaded {len(discoveries)} discoveries")
    return discoveries

def create_consolidated_dataset(discoveries):
    """Create consolidated dataset with summary"""
    
    # Calculate statistics
    scores = [d.get('combined_score', d.get('score', 0)) for d in discoveries if d.get('combined_score') or d.get('score')]
    avg_score = sum(scores) / len(scores) if scores else 0
    max_score = max(scores) if scores else 0
    min_score = min(scores) if scores else 0
    
    # Get unique SMILES
    unique_smiles = set()
    for d in discoveries:
        smiles = d.get('smiles', '')
        if smiles:
            unique_smiles.add(smiles)
    
    # Create consolidated dataset
    dataset = {
        "discovery_summary": {
            "total_discoveries": len(discoveries),
            "generated_at": datetime.now().isoformat(),
            "discovery_campaign": "Continuous Autonomous Discovery - Latest Update",
            "statistics": {
                "avg_score": avg_score,
                "max_score": max_score,
                "min_score": min_score,
                "unique_molecules": len(unique_smiles),
                "score_field_mapped": True
            },
            "data_version": "3.0",
            "latest_update": True,
            "individual_files_consolidated": len(glob('continuous_chemistry_discoveries/discoveries/*.json'))
        },
        "discoveries": []
    }
    
    # Process discoveries and ensure score mapping
    for discovery in discoveries:
        # Ensure score field exists for Streamlit compatibility
        if 'combined_score' in discovery and 'score' not in discovery:
            discovery['score'] = discovery['combined_score']
        elif 'score' not in discovery and 'combined_score' not in discovery:
            discovery['score'] = 0.0
            
        dataset["discoveries"].append(discovery)
    
    return dataset

def update_data_files(dataset):
    """Update all the data files that Streamlit apps expect"""
    
    files_to_update = [
        "results/chemistry_discoveries.json",
        "cloud_data_snapshot.json", 
        "results/overnight_discovery_mega_dataset.json"
    ]
    
    for file_path in files_to_update:
        print(f"📝 Updating {file_path}...")
        try:
            with open(file_path, 'w') as f:
                json.dump(dataset, f, indent=2)
            print(f"   ✅ Updated {file_path} with {len(dataset['discoveries'])} discoveries")
        except Exception as e:
            print(f"   ❌ Failed to update {file_path}: {e}")

def update_problem_solution_data():
    """Update problem-solution data if analyzer exists"""
    try:
        print("🎯 Updating problem-solution analysis...")
        os.system("python3 problem_solution_analyzer.py > /dev/null 2>&1")
        
        # Copy to cloud files if analysis was successful
        if os.path.exists("problem_solution_analysis/complete_analysis.json"):
            os.system("cp problem_solution_analysis/complete_analysis.json cloud_problem_solution_data.json")
            os.system("cp problem_solution_analysis/complete_analysis.json results/problem_solution_summary.json")
            print("   ✅ Problem-solution analysis updated")
        else:
            print("   ⚠️ Problem-solution analysis not available")
    except Exception as e:
        print(f"   ⚠️ Problem-solution update failed: {e}")

if __name__ == "__main__":
    print("🚀 UPDATING DISCOVERY DATA WITH LATEST FINDINGS...")
    print()
    
    # Load all individual discoveries
    discoveries = load_individual_discoveries()
    
    if not discoveries:
        print("❌ No discoveries found to consolidate")
        sys.exit(1)
    
    # Create consolidated dataset
    print("📊 Creating consolidated dataset...")
    dataset = create_consolidated_dataset(discoveries)
    
    print(f"✅ Dataset created:")
    print(f"   🧬 Total discoveries: {dataset['discovery_summary']['total_discoveries']}")
    print(f"   📊 Average score: {dataset['discovery_summary']['statistics']['avg_score']:.3f}")
    print(f"   🏆 Max score: {dataset['discovery_summary']['statistics']['max_score']:.3f}")
    print(f"   🔬 Unique molecules: {dataset['discovery_summary']['statistics']['unique_molecules']}")
    print()
    
    # Update all data files
    print("📦 Updating data files...")
    update_data_files(dataset)
    print()
    
    # Update problem-solution analysis
    update_problem_solution_data()
    print()
    
    print("🎉 DATA UPDATE COMPLETE!")
    print(f"   📊 {len(discoveries)} discoveries now available in all Streamlit apps")
    print(f"   🌐 Cloud deployment ready with latest data")
    print(f"   🔄 All file names preserved for compatibility")
