#!/usr/bin/env python3
"""
Fix for molecular generation to eliminate duplicates and ensure quality molecules
"""

import sys
sys.path.append('agents/alchemist')
sys.path.append('core')

from rdkit import Chem
from rdkit.Chem import Descriptors
import json
import os
from typing import Set, List, Dict, Any

def validate_molecule_quality(smiles: str) -> tuple[bool, str]:
    """Validate if a molecule meets quality criteria"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES"
        
        # Basic quality checks
        num_atoms = mol.GetNumAtoms()
        num_heavy_atoms = mol.GetNumHeavyAtoms()
        
        # Minimum complexity requirements
        if num_atoms < 3:
            return False, "Too simple (< 3 atoms)"
        
        if num_heavy_atoms < 2:
            return False, "Too simple (< 2 heavy atoms)"
        
        # Maximum complexity (to avoid huge molecules)
        if num_atoms > 100:
            return False, "Too complex (> 100 atoms)"
        
        # Check for reasonable molecular weight
        mw = Descriptors.MolWt(mol)
        if mw < 30:  # Lighter than formaldehyde
            return False, f"Too light (MW: {mw:.1f})"
        
        if mw > 1000:  # Too heavy for drug-like
            return False, f"Too heavy (MW: {mw:.1f})"
        
        # Check for valid connectivity
        if mol.GetNumBonds() == 0 and num_atoms > 1:
            return False, "No bonds between atoms"
        
        # Canonical SMILES check (eliminates tautomers and stereoisomers)
        canonical = Chem.MolToSmiles(mol)
        if canonical != smiles:
            return False, f"Non-canonical SMILES (canonical: {canonical})"
        
        return True, "Valid"
        
    except Exception as e:
        return False, f"Error: {e}"

def clean_discovery_files():
    """Clean up discovery files by removing duplicates and invalid molecules"""
    discovery_dir = "continuous_chemistry_discoveries/discoveries"
    if not os.path.exists(discovery_dir):
        print(f"❌ Discovery directory not found: {discovery_dir}")
        return
    
    files = [f for f in os.listdir(discovery_dir) if f.endswith('.json')]
    print(f"🔍 Found {len(files)} discovery files")
    
    all_discoveries = []
    seen_smiles: Set[str] = set()
    valid_count = 0
    invalid_count = 0
    duplicate_count = 0
    
    for filename in files:
        filepath = os.path.join(discovery_dir, filename)
        try:
            with open(filepath, 'r') as f:
                discovery = json.load(f)
            
            smiles = discovery.get('smiles', '')
            if not smiles:
                invalid_count += 1
                continue
            
            # Check for duplicates
            if smiles in seen_smiles:
                duplicate_count += 1
                print(f"🗑️ Removing duplicate: {smiles}")
                os.remove(filepath)  # Delete duplicate file
                continue
            
            # Validate molecule quality
            is_valid, reason = validate_molecule_quality(smiles)
            if not is_valid:
                invalid_count += 1
                print(f"🗑️ Removing invalid molecule ({reason}): {smiles}")
                os.remove(filepath)  # Delete invalid file
                continue
            
            # Keep valid, unique molecule
            seen_smiles.add(smiles)
            all_discoveries.append(discovery)
            valid_count += 1
            
        except Exception as e:
            print(f"❌ Error processing {filename}: {e}")
            invalid_count += 1
    
    print(f"\n📊 Cleanup Results:")
    print(f"✅ Valid unique molecules: {valid_count}")
    print(f"🗑️ Duplicates removed: {duplicate_count}")
    print(f"❌ Invalid molecules removed: {invalid_count}")
    print(f"📋 Total processed: {len(files)}")
    
    # Show some examples of remaining molecules
    if all_discoveries:
        print(f"\n🧬 Sample of remaining valid molecules:")
        sorted_discoveries = sorted(all_discoveries, key=lambda x: x.get('combined_score', 0), reverse=True)
        for i, d in enumerate(sorted_discoveries[:10]):
            smiles = d['smiles']
            score = d.get('combined_score', 0)
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mw = Descriptors.MolWt(mol)
                print(f"  {i+1}. {smiles} (Score: {score:.3f}, MW: {mw:.1f})")
    
    return all_discoveries

def analyze_discovery_diversity():
    """Analyze the diversity of discovered molecules"""
    discovery_dir = "continuous_chemistry_discoveries/discoveries"
    if not os.path.exists(discovery_dir):
        return
    
    files = [f for f in os.listdir(discovery_dir) if f.endswith('.json')]
    
    molecular_formulas = {}
    smiles_set = set()
    scores = []
    
    for filename in files:
        filepath = os.path.join(discovery_dir, filename)
        try:
            with open(filepath, 'r') as f:
                discovery = json.load(f)
            
            smiles = discovery.get('smiles', '')
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Get molecular formula
                    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                    molecular_formulas[formula] = molecular_formulas.get(formula, 0) + 1
                    smiles_set.add(smiles)
                    scores.append(discovery.get('combined_score', 0))
        except:
            continue
    
    print(f"\n🔬 Discovery Diversity Analysis:")
    print(f"📊 Unique SMILES: {len(smiles_set)}")
    print(f"🧪 Unique molecular formulas: {len(molecular_formulas)}")
    if scores:
        import statistics
        print(f"📈 Score range: {min(scores):.3f} - {max(scores):.3f}")
        print(f"📊 Average score: {statistics.mean(scores):.3f}")
    
    print(f"\n🏆 Most common molecular formulas:")
    sorted_formulas = sorted(molecular_formulas.items(), key=lambda x: x[1], reverse=True)
    for formula, count in sorted_formulas[:10]:
        print(f"  {formula}: {count} molecules")

if __name__ == "__main__":
    print("🧬 FoTChemistry Molecular Quality Control")
    print("=" * 50)
    
    # Clean up discovery files
    valid_discoveries = clean_discovery_files()
    
    # Analyze diversity
    analyze_discovery_diversity()
    
    # Update the Streamlit export
    if valid_discoveries:
        print(f"\n📤 Updating Streamlit export...")
        
        # Sort by score
        valid_discoveries.sort(key=lambda x: x.get('combined_score', 0), reverse=True)
        
        # Create export data
        export_data = {
            'export_timestamp': '2025-09-28T10:54:00Z',
            'total_discoveries': len(valid_discoveries),
            'discoveries': valid_discoveries,
            'recent_molecules': valid_discoveries[:20],
            'statistics': {
                'total_molecules': len(valid_discoveries),
                'avg_score': sum(d.get('combined_score', 0) for d in valid_discoveries) / len(valid_discoveries) if valid_discoveries else 0,
                'max_score': max(d.get('combined_score', 0) for d in valid_discoveries) if valid_discoveries else 0,
                'total_reactions': 0,
                'total_measurements': 0,
                'active_claims': 0
            }
        }
        
        # Transform for Streamlit format
        streamlit_discoveries = []
        for d in valid_discoveries:
            drug_likeness = d.get('drug_likeness', {})
            if isinstance(drug_likeness, dict):
                drug_likeness_display = drug_likeness
            else:
                drug_likeness_display = {'passes_lipinski': False}
                
            streamlit_discovery = {
                'id': d.get('discovery_id', 'unknown'),
                'smiles': d.get('smiles', ''),
                'score': d.get('combined_score', 0.0),
                'drug_likeness': drug_likeness_display,
                'safety_score': d.get('safety_assessment', {}).get('safety_score', 1.0),
                'quantum_coherence': d.get('quantum_measurements', {}).get('coherence', 0.0),
                'timestamp': d.get('discovery_timestamp', ''),
                'properties': {
                    'molecular_weight': d.get('molecular_properties', {}).get('molecular_weight', 0.0),
                    'logp': d.get('molecular_properties', {}).get('logp', 0.0),
                    'molecular_formula': d.get('molecular_properties', {}).get('molecular_formula', 'N/A')
                }
            }
            streamlit_discoveries.append(streamlit_discovery)
        
        export_data['discoveries'] = streamlit_discoveries
        export_data['recent_molecules'] = streamlit_discoveries[:20]
        
        # Save to results
        os.makedirs('results', exist_ok=True)
        with open('results/chemistry_discoveries.json', 'w') as f:
            json.dump(export_data, f, indent=2)
        
        print(f"✅ Updated Streamlit export with {len(valid_discoveries)} clean discoveries")
        print(f"📊 Average score: {export_data['statistics']['avg_score']:.3f}")
        print(f"🏆 Max score: {export_data['statistics']['max_score']:.3f}")
    
    print(f"\n🎉 Molecular quality control complete!")
    print(f"🌐 View results at: http://localhost:8505")
