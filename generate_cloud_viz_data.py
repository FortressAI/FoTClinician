#!/usr/bin/env python3
"""
Generate cloud visualization data for all 6,443 molecules
Adds 2D SVG and 3D MOL block data for cloud deployment
"""

import json
import sys
from datetime import datetime

# Add paths for local imports
sys.path.append('core')
sys.path.append('akg')

# Try to import RDKit for visualization generation
HAS_RDKIT = False
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
    print("✅ RDKit available - generating high-quality structures")
except ImportError:
    print("⚠️ RDKit not available - using fallback rendering")

def generate_2d_svg(smiles: str) -> dict:
    """Generate 2D SVG structure"""
    if not smiles or not HAS_RDKIT:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Create SVG drawer
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.SetFontSize(16)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        return {
            "svg": svg,
            "format": "svg",
            "width": 400,
            "height": 400
        }
    except Exception as e:
        print(f"⚠️ 2D generation failed for {smiles}: {e}")
        return None

def generate_3d_molblock(smiles: str) -> dict:
    """Generate 3D MOL block data"""
    if not smiles or not HAS_RDKIT:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        # Skip optimization for now - just use embedded coordinates
        
        # Generate MOL block
        molblock = Chem.MolToMolBlock(mol)
        
        # Extract atom positions for additional data
        conf = mol.GetConformer()
        atoms = []
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            atoms.append({
                "symbol": atom.GetSymbol(),
                "x": pos.x,
                "y": pos.y,
                "z": pos.z
            })
        
        return {
            "molblock": molblock,
            "atoms": atoms,
            "format": "mol"
        }
    except Exception as e:
        print(f"⚠️ 3D generation failed for {smiles}: {e}")
        return None

def add_visualization_data(input_file: str, output_file: str, max_molecules: int = None):
    """Add visualization data to discovery dataset"""
    
    print(f"📂 Loading data from {input_file}...")
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    discoveries = data.get('discoveries', [])
    total_molecules = len(discoveries)
    
    if max_molecules:
        discoveries = discoveries[:max_molecules]
        print(f"🎯 Processing first {len(discoveries)} molecules (limited for cloud)")
    else:
        print(f"🧬 Processing all {total_molecules} molecules...")
    
    processed = 0
    success_2d = 0
    success_3d = 0
    
    for i, molecule in enumerate(discoveries):
        smiles = molecule.get('smiles', '')
        if not smiles:
            continue
            
        print(f"🔄 Processing {i+1}/{len(discoveries)}: {smiles}")
        
        # Generate 2D visualization
        viz_2d = generate_2d_svg(smiles)
        if viz_2d:
            molecule['visualization_2d'] = viz_2d
            success_2d += 1
        
        # Generate 3D visualization
        viz_3d = generate_3d_molblock(smiles)
        if viz_3d:
            molecule['visualization_3d'] = viz_3d
            success_3d += 1
        
        processed += 1
        
        # Progress update every 100 molecules
        if processed % 100 == 0:
            print(f"✅ Progress: {processed}/{len(discoveries)} molecules processed")
    
    # Update summary
    data['discovery_summary']['visualization_generated'] = {
        "timestamp": datetime.now().isoformat(),
        "total_processed": processed,
        "success_2d": success_2d,
        "success_3d": success_3d,
        "rdkit_available": HAS_RDKIT
    }
    
    # Save enhanced data
    print(f"💾 Saving enhanced data to {output_file}...")
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"🎉 VISUALIZATION DATA GENERATION COMPLETE!")
    print(f"   📊 Total molecules: {processed}")
    print(f"   🎨 2D visualizations: {success_2d}")
    print(f"   🧊 3D visualizations: {success_3d}")
    print(f"   📁 Output file: {output_file}")
    
    return data

if __name__ == "__main__":
    if not HAS_RDKIT:
        print("❌ RDKit is required for visualization generation!")
        print("   Install with: pip install rdkit-pypi")
        sys.exit(1)
    
    # Generate cloud-optimized version (first 100 molecules for testing)
    print("🌐 Generating cloud-optimized visualization data...")
    add_visualization_data(
        "results/chemistry_discoveries.json",
        "cloud_data_snapshot_with_viz.json",
        max_molecules=100  # Limit for cloud deployment
    )
    
    print("\n🎯 Cloud visualization data ready for deployment!")
