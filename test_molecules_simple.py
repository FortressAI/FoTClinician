#!/usr/bin/env python3
"""
Simple test to prove 2D and 3D molecular visualization works
"""

import streamlit as st
import streamlit.components.v1 as components

# Try imports
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    import stmol
    HAS_STMOL = True
except ImportError:
    HAS_STMOL = False

try:
    import py3Dmol
    HAS_PY3DMOL = True
except ImportError:
    HAS_PY3DMOL = False

def test_2d_molecule(smiles):
    """Test 2D molecule rendering"""
    st.subheader(f"2D Structure: {smiles}")
    
    if not HAS_RDKIT:
        st.error("RDKit not available")
        return
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.error(f"Invalid SMILES: {smiles}")
            return
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Create SVG
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.SetFontSize(16)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        # Clean SVG (remove XML declaration)
        if svg.startswith('<?xml'):
            svg = svg[svg.find('<svg'):]
        
        st.markdown(svg, unsafe_allow_html=True)
        st.success("✅ 2D rendering successful")
        
    except Exception as e:
        st.error(f"2D rendering failed: {e}")

def test_3d_molecule(smiles):
    """Test 3D molecule rendering"""
    st.subheader(f"3D Structure: {smiles}")
    
    if not HAS_RDKIT:
        st.error("RDKit not available")
        return
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.error(f"Invalid SMILES: {smiles}")
            return
        
        # Generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)
        
        # Get MOL block
        molblock = Chem.MolToMolBlock(mol)
        
        # Try stmol first (CORRECTED VERSION)
        if HAS_STMOL and HAS_PY3DMOL:
            try:
                st.write("Trying stmol with py3Dmol...")
                # CORRECT WAY: Create py3Dmol viewer first
                viewer = py3Dmol.view(width=400, height=400)
                viewer.addModel(molblock, 'mol')
                viewer.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.3}})
                viewer.setBackgroundColor('white')
                viewer.zoomTo()
                
                # Pass the VIEWER object to stmol, not the molblock string
                stmol.showmol(viewer, height=400, width=400)
                st.success("✅ 3D rendering with stmol successful")
                return
                
            except Exception as e:
                st.warning(f"stmol failed: {e}")
        
        # Fallback to direct py3Dmol
        if HAS_PY3DMOL:
            try:
                st.write("Trying direct py3Dmol...")
                viewer = py3Dmol.view(width=400, height=400)
                viewer.addModel(molblock, 'mol')
                viewer.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.3}})
                viewer.setBackgroundColor('white')
                viewer.zoomTo()
                
                # Get HTML and display
                html_content = viewer._make_html()
                components.html(html_content, height=400, width=400)
                st.success("✅ 3D rendering with py3Dmol successful")
                return
                
            except Exception as e:
                st.warning(f"py3Dmol failed: {e}")
        
        # Final fallback
        st.info("3D libraries not available - showing MOL block")
        st.code(molblock, language='text')
        
    except Exception as e:
        st.error(f"3D rendering failed: {e}")

def main():
    st.title("🧪 Simple Molecular Visualization Test")
    
    # Show what's available
    st.subheader("Available Libraries")
    st.write(f"RDKit: {'✅' if HAS_RDKIT else '❌'}")
    st.write(f"stmol: {'✅' if HAS_STMOL else '❌'}")
    st.write(f"py3Dmol: {'✅' if HAS_PY3DMOL else '❌'}")
    
    # Test molecules
    test_molecules = [
        "CCO",      # Ethanol
        "CCCO",     # Propanol
        "CC(=O)O",  # Acetic acid
        "C1=CC=CC=C1"  # Benzene
    ]
    
    selected_smiles = st.selectbox("Choose a molecule to test:", test_molecules)
    
    col1, col2 = st.columns(2)
    
    with col1:
        test_2d_molecule(selected_smiles)
    
    with col2:
        test_3d_molecule(selected_smiles)

if __name__ == "__main__":
    main()
