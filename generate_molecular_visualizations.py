#!/usr/bin/env python3
"""
Generate molecular visualizations for cloud deployment
Pre-renders all 2D/3D structures and saves them as data files
"""

import json
import os
import sys
import base64
from pathlib import Path

# Add paths for local imports
sys.path.append('core')
sys.path.append('akg')

# Try to import RDKit for high-quality rendering
HAS_RDKIT = False
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
    print("✅ RDKit available - generating high-quality structures")
except ImportError:
    print("⚠️ RDKit not available - using fallback rendering")

# Import our custom visualization
from simple_mol_viz import smiles_to_simple_svg, create_3d_fallback_html

def generate_2d_svg(smiles: str) -> str:
    """Generate 2D SVG structure"""
    if not smiles:
        return None
    
    # Try RDKit first
    if HAS_RDKIT:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
                drawer.SetFontSize(16)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                return svg
        except Exception as e:
            print(f"RDKit failed for {smiles}: {e}")
    
    # Fallback to custom SVG
    try:
        svg = smiles_to_simple_svg(smiles, 400, 400)
        return svg
    except Exception as e:
        print(f"Custom SVG failed for {smiles}: {e}")
        return None

def generate_3d_data(smiles: str) -> dict:
    """Generate 3D molecular data"""
    if not smiles:
        return None
    
    # Try RDKit for real 3D coordinates
    if HAS_RDKIT:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.UFFOptimizeMolecule(mol)
                
                # Get MOL block for 3D viewers
                molblock = Chem.MolToMolBlock(mol)
                
                # Also get XYZ coordinates
                conf = mol.GetConformer()
                atoms = []
                for i, atom in enumerate(mol.GetAtoms()):
                    pos = conf.GetAtomPosition(i)
                    atoms.append({
                        'symbol': atom.GetSymbol(),
                        'x': pos.x,
                        'y': pos.y,
                        'z': pos.z
                    })
                
                return {
                    'molblock': molblock,
                    'atoms': atoms,
                    'has_3d': True,
                    'method': 'rdkit'
                }
        except Exception as e:
            print(f"RDKit 3D failed for {smiles}: {e}")
    
    # Fallback to HTML representation
    try:
        html = create_3d_fallback_html(smiles)
        return {
            'html': html,
            'has_3d': False,
            'method': 'fallback'
        }
    except Exception as e:
        print(f"3D fallback failed for {smiles}: {e}")
        return None

def process_discovery_data():
    """Process all discovery data and add visualizations"""
    
    # Load discovery data
    discovery_files = [
        'cloud_data_snapshot.json',
        'results/chemistry_discoveries.json'
    ]
    
    for discovery_file in discovery_files:
        if not os.path.exists(discovery_file):
            print(f"⚠️ Discovery file not found: {discovery_file}")
            continue
            
        print(f"🔄 Processing {discovery_file}...")
        
        with open(discovery_file, 'r') as f:
            data = json.load(f)
        
        # Get discoveries list
        discoveries = data.get('discoveries', [])
        if not discoveries:
            print(f"❌ No discoveries found in {discovery_file}")
            continue
        
        print(f"📊 Found {len(discoveries)} molecules to process")
        
        # Process each molecule
        processed_count = 0
        for i, discovery in enumerate(discoveries):
            smiles = discovery.get('smiles', '')
            if not smiles:
                continue
            
            print(f"🧬 Processing {i+1}/{len(discoveries)}: {smiles}")
            
            # Generate 2D structure
            svg_2d = generate_2d_svg(smiles)
            if svg_2d:
                discovery['visualization_2d'] = {
                    'svg': svg_2d,
                    'format': 'svg',
                    'width': 400,
                    'height': 400
                }
            
            # Generate 3D structure
            data_3d = generate_3d_data(smiles)
            if data_3d:
                discovery['visualization_3d'] = data_3d
            
            processed_count += 1
        
        # Save enhanced data
        output_file = discovery_file.replace('.json', '_with_viz.json')
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"✅ Saved {processed_count} molecules with visualizations to {output_file}")
        
        # Create a summary
        viz_summary = {
            'total_molecules': len(discoveries),
            'with_2d': len([d for d in discoveries if 'visualization_2d' in d]),
            'with_3d': len([d for d in discoveries if 'visualization_3d' in d]),
            'rdkit_available': HAS_RDKIT,
            'source_file': discovery_file,
            'output_file': output_file
        }
        
        summary_file = output_file.replace('.json', '_summary.json')
        with open(summary_file, 'w') as f:
            json.dump(viz_summary, f, indent=2)
        
        print(f"📋 Summary saved to {summary_file}")

def create_cloud_optimized_requirements():
    """Create requirements file without problematic packages"""
    
    cloud_requirements = [
        "streamlit>=1.28.0",
        "pandas>=1.5.0", 
        "numpy>=1.24.0",
        "plotly>=5.15.0",
        "requests>=2.31.0"
    ]
    
    with open('requirements_cloud_optimized.txt', 'w') as f:
        f.write('\n'.join(cloud_requirements))
    
    print("✅ Created requirements_cloud_optimized.txt (no chemistry packages)")

if __name__ == "__main__":
    print("🧬 FoTChemistry Molecular Visualization Generator")
    print("=" * 60)
    
    # Generate visualizations
    process_discovery_data()
    
    # Create optimized requirements
    create_cloud_optimized_requirements()
    
    print("\n🎉 Visualization generation complete!")
    print("\n📋 Next steps for cloud deployment:")
    print("1. Use requirements_cloud_optimized.txt as requirements.txt")
    print("2. Update Streamlit app to use pre-generated visualizations")
    print("3. Commit visualization data files to repository")
    print("4. Deploy to Streamlit Cloud")
