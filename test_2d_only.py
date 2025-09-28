#!/usr/bin/env python3
"""
Simple test for 2D molecular visualization only
"""

import streamlit as st

try:
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

def render_2d_molecule(smiles):
    """Simple 2D molecule rendering"""
    if not HAS_RDKIT:
        st.error("RDKit not available")
        return
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.error(f"Invalid SMILES: {smiles}")
            return
        
        st.success(f"✅ Molecule created from SMILES: {smiles}")
        
        # Generate 2D coordinates
        rdDepictor.Compute2DCoords(mol)
        st.success("✅ 2D coordinates computed")
        
        # Create SVG
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.SetFontSize(16)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        st.success(f"✅ SVG generated, length: {len(svg)}")
        st.write(f"Has bonds: {'bond-' in svg}")
        st.write(f"Has XML declaration: {svg.startswith('<?xml')}")
        
        # Clean SVG
        cleaned_svg = svg
        if cleaned_svg.startswith('<?xml'):
            cleaned_svg = cleaned_svg[cleaned_svg.find('<svg'):]
            st.success(f"✅ SVG cleaned, new length: {len(cleaned_svg)}")
        
        # Show raw SVG first
        st.subheader("Raw SVG (first 500 chars):")
        st.code(cleaned_svg[:500], language='xml')
        
        # Try displaying it
        st.subheader("Rendered 2D Structure:")
        st.markdown(cleaned_svg, unsafe_allow_html=True)
        
        # Alternative: try with components
        st.subheader("Alternative: Using st.components.v1.html:")
        import streamlit.components.v1 as components
        components.html(cleaned_svg, height=400, width=400)
        
    except Exception as e:
        st.error(f"Error: {e}")
        import traceback
        st.code(traceback.format_exc())

def main():
    st.title("🧪 2D Molecular Visualization Debug")
    
    st.write(f"RDKit available: {'✅' if HAS_RDKIT else '❌'}")
    
    # Test molecules
    test_molecules = ["CCO", "CCCO", "CC(=O)O", "C1=CC=CC=C1"]
    selected = st.selectbox("Choose molecule:", test_molecules)
    
    if st.button("Render 2D Structure"):
        render_2d_molecule(selected)

if __name__ == "__main__":
    main()
