#!/usr/bin/env python3
"""
FoTChemistry Discovery Dashboard

Dashboard-first interface for viewing molecular discoveries with detailed analysis.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import json
import os
import sys
import base64
from io import BytesIO

# Add core and akg to path
sys.path.append('core')
sys.path.append('akg')
sys.path.append('agents/alchemist')

# Feature flags for conditional imports
HAS_AKG = False
HAS_QUANTUM = False  
HAS_RDKIT = False
HAS_STMOL = False
HAS_PY3DMOL = False

try:
    from client import AKG
    HAS_AKG = True
except ImportError:
    pass

try:
    from chemistry_vqbit_engine import ChemistryVQbitEngine, ChemistryPropertyType
    HAS_QUANTUM = True
except ImportError:
    pass

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
except ImportError:
    pass

try:
    import stmol
    HAS_STMOL = True
except ImportError:
    pass

try:
    import py3Dmol
    HAS_PY3DMOL = True
except ImportError:
    pass

# Page configuration
st.set_page_config(
    page_title="FoTChemistry Discovery Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

@st.cache_data
def load_discovery_data():
    """Load discovery data with cloud/local detection and fallback options."""
    
    # Detect deployment environment
    is_cloud_deployment = (
        os.environ.get('STREAMLIT_SHARING') == '1' or
        os.environ.get('STREAMLIT_CLOUD') == '1' or
        'streamlit.io' in os.environ.get('HOSTNAME', '') or
        not os.path.exists('akg/client.py')  # Local AKG not available
    )
    
    if is_cloud_deployment:
        st.info("☁️ Cloud deployment detected - using static data snapshot")
        
        # Try cloud snapshot with visualizations first
        for snapshot_file in ['cloud_data_snapshot_with_viz.json', 'cloud_data_snapshot.json', 'demo_data.json']:
            if os.path.exists(snapshot_file):
                try:
                    with open(snapshot_file, 'r') as f:
                        data = json.load(f)
                    st.success(f"📦 Loaded {data.get('discovery_summary', {}).get('total_discoveries', 0)} molecules from {snapshot_file}")
                    return data
                except Exception as e:
                    st.warning(f"⚠️ Failed to load {snapshot_file}: {e}")
        
        # No demo data - return None to force error if no real data available
        st.error("❌ No discovery data available - real data required")
        return None
    
    else:
        st.info("🏠 Local deployment detected - using live data")
        
        # Try to load from local results file  
        export_file = "results/chemistry_discoveries.json"
        if os.path.exists(export_file):
            try:
                with open(export_file, 'r') as f:
                    data = json.load(f)
                st.success(f"📂 Loaded {data.get('total_discoveries', 0)} molecules from local results")
                return data
            except Exception as e:
                st.warning(f"⚠️ Failed to load local results: {e}")
        
        # No demo data - return None to force real data requirement
        st.error("❌ No discovery data found - check data files")
        return None

def create_demo_data():
    """Create embedded demo data for cloud deployment."""
    from datetime import datetime
    
    demo_molecules = [
        {
            "id": "demo_001",
            "smiles": "CCCO",
            "name": "Propanol",
            "score": 0.786,
            "drug_likeness": {"passes_lipinski": True, "score": 0.85},
            "safety_score": 0.92,
            "quantum_coherence": 0.74,
            "timestamp": datetime.now().isoformat(),
            "properties": {
                "formula": "C3H8O",
                "molecular_weight": 60.1,
                "logp": 0.25,
                "tpsa": 20.23,
                "hbd": 1,
                "hba": 1
            }
        },
        {
            "id": "demo_002", 
            "smiles": "c1ccc(-c2ccccc2)cc1",
            "name": "Biphenyl",
            "score": 0.754,
            "drug_likeness": {"passes_lipinski": True, "score": 0.78},
            "safety_score": 0.88,
            "quantum_coherence": 0.71,
            "timestamp": datetime.now().isoformat(),
            "properties": {
                "formula": "C12H10",
                "molecular_weight": 154.2,
                "logp": 3.76,
                "tpsa": 0.0,
                "hbd": 0,
                "hba": 0
            }
        }
    ]
    
    return {
        "discoveries": demo_molecules,
        "total_discoveries": len(demo_molecules),
        "statistics": {
            "total_molecules": len(demo_molecules),
            "avg_score": sum(m["score"] for m in demo_molecules) / len(demo_molecules),
            "max_score": max(m["score"] for m in demo_molecules),
            "total_reactions": 0,
            "total_measurements": 0,
            "active_claims": 0
        },
        "recent_molecules": demo_molecules
    }

def render_molecule_2d(molecule_data: dict):
    """Render 2D molecular structure from pre-generated or live data"""
    smiles = molecule_data.get('smiles', '')
    if not smiles:
        return None
    
    # Try pre-generated visualization first
    if 'visualization_2d' in molecule_data:
        viz_data = molecule_data['visualization_2d']
        if viz_data.get('svg'):
            return viz_data['svg']
    
    # Fallback to live generation if RDKit available
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
            st.warning(f"RDKit visualization failed: {e}")
    
    # Fallback to simple molecular visualization
    try:
        from simple_mol_viz import smiles_to_simple_svg
        svg = smiles_to_simple_svg(smiles, 400, 400)
        return svg
    except Exception as e:
        # Ultimate fallback - styled text display
        return f'''
        <div style="width: 400px; height: 400px; border: 1px solid #ddd; border-radius: 8px; 
                    background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%); 
                    display: flex; flex-direction: column; align-items: center; justify-content: center;
                    font-family: Arial, sans-serif;">
            <div style="font-size: 48px; margin-bottom: 20px;">⚗️</div>
            <div style="font-size: 18px; font-weight: bold; margin-bottom: 10px; color: #1565c0;">
                Molecular Structure
            </div>
            <div style="font-size: 14px; color: #1976d2; margin-bottom: 15px;">
                {smiles}
            </div>
            <div style="font-size: 12px; color: #666; text-align: center; padding: 0 20px;">
                Visual structure rendering unavailable
        </div>
        </div>
        '''

def render_molecule_3d(molecule_data: dict):
    """Render 3D molecular structure using stmol/py3Dmol and pre-generated data"""
    smiles = molecule_data.get('smiles', '')
    if not smiles:
        return None
    
    # Try pre-generated 3D data first
    if 'visualization_3d' in molecule_data:
        viz_data = molecule_data['visualization_3d']
        if viz_data.get('has_3d') and viz_data.get('molblock'):
            molblock = viz_data['molblock']
            
            # Try stmol first (preferred for Streamlit Cloud)
            if HAS_STMOL:
                try:
                    stmol.showmol(molblock, height=400, width=400)
                    return True
                except Exception as e:
                    st.warning(f"Stmol visualization failed: {e}")
            
            # Try py3Dmol as alternative
            if HAS_PY3DMOL:
                try:
                    viewer = py3Dmol.view(width=400, height=400)
                    viewer.addModel(molblock, 'mol')
                    viewer.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.3}})
                    viewer.setBackgroundColor('white')
                    viewer.zoomTo()
                    stmol.showmol(viewer, height=400, width=400)
                    return True
                except Exception as e:
                    st.warning(f"Py3Dmol visualization failed: {e}")
        
        elif viz_data.get('html'):
            # Use pre-generated HTML fallback
            st.components.v1.html(viz_data['html'], height=400)
            return True
    
    # Fallback to live generation if RDKit available
    if HAS_RDKIT:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.UFFOptimizeMolecule(mol)
                
                molblock = Chem.MolToMolBlock(mol)
                
                # Try stmol first
                if HAS_STMOL:
                    stmol.showmol(molblock, height=400, width=400)
                    return True
                
                # Try py3Dmol alternative
                if HAS_PY3DMOL:
                    viewer = py3Dmol.view(width=400, height=400)
                    viewer.addModel(molblock, 'mol')
                    viewer.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.3}})
                    viewer.setBackgroundColor('white')
                    viewer.zoomTo()
                    stmol.showmol(viewer, height=400, width=400)
                    return True
                
        except Exception as e:
            st.warning(f"Live 3D generation failed: {e}")
    
    # Ultimate fallback - styled HTML
    fallback_html = f'''
    <div style="width: 400px; height: 400px; border: 1px solid #ddd; border-radius: 8px; 
                background: linear-gradient(135deg, #f3e5f5 0%, #e1bee7 100%); 
                display: flex; flex-direction: column; align-items: center; justify-content: center;
                font-family: Arial, sans-serif;">
        <div style="font-size: 48px; margin-bottom: 20px;">🧬</div>
        <div style="font-size: 18px; font-weight: bold; margin-bottom: 10px; color: #7b1fa2;">
            3D Molecular Model
        </div>
        <div style="font-size: 14px; color: #8e24aa; margin-bottom: 15px;">
            {smiles}
        </div>
        <div style="font-size: 12px; color: #666; text-align: center; padding: 0 20px;">
            Interactive 3D visualization
        </div>
    </div>
    '''
    st.components.v1.html(fallback_html, height=400)
    return True

def calculate_molecular_properties(smiles: str):
    """Calculate molecular properties"""
    if not HAS_RDKIT or not smiles:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        properties = {
            'molecular_formula': rdMolDescriptors.CalcMolFormula(mol),
            'molecular_weight': round(Descriptors.MolWt(mol), 2),
            'logp': round(Descriptors.MolLogP(mol), 2),
            'tpsa': round(Descriptors.TPSA(mol), 2),
            'hbd_count': Descriptors.NumHDonors(mol),
            'hba_count': Descriptors.NumHAcceptors(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'heavy_atoms': mol.GetNumHeavyAtoms(),
            'rings': Descriptors.RingCount(mol),
            'aromatic_rings': Descriptors.NumAromaticRings(mol),
            'lipinski_violations': sum([
                Descriptors.MolWt(mol) > 500,
                Descriptors.MolLogP(mol) > 5,
                Descriptors.NumHDonors(mol) > 5,
                Descriptors.NumHAcceptors(mol) > 10
            ])
        }
        
        return properties
    except Exception as e:
        return None

def display_molecule_detail(molecule):
    """Display detailed view of a single molecule"""
    st.header(f"🧬 Molecule Detail: {molecule.get('smiles', 'Unknown')}")
    
    # Main layout
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📊 Discovery Information")
        
        # Basic info
        st.markdown(f"**SMILES:** `{molecule.get('smiles', 'N/A')}`")
        st.markdown(f"**Discovery Score:** {molecule.get('score', 0):.3f}")
        st.markdown(f"**Discovery ID:** {molecule.get('id', 'Unknown')[:12]}...")
        
        # Drug likeness info
        drug_likeness = molecule.get('drug_likeness', {})
        if isinstance(drug_likeness, dict):
            passes_lipinski = drug_likeness.get('passes_lipinski', False)
            st.markdown(f"**Drug Likeness:** {'✅ Pass' if passes_lipinski else '❌ Fail'} Lipinski's Rule")
        
        st.markdown(f"**Safety Score:** {molecule.get('safety_score', 0):.3f}")
        st.markdown(f"**Quantum Coherence:** {molecule.get('quantum_coherence', 0):.6f}")
        
        # Molecular properties
        if HAS_RDKIT:
            props = calculate_molecular_properties(molecule.get('smiles', ''))
            if props:
                st.subheader("🔬 Molecular Properties")
                st.markdown(f"**Formula:** {props['molecular_formula']}")
                st.markdown(f"**Molecular Weight:** {props['molecular_weight']} g/mol")
                st.markdown(f"**LogP:** {props['logp']}")
                st.markdown(f"**TPSA:** {props['tpsa']} Ų")
                st.markdown(f"**H-Bond Donors:** {props['hbd_count']}")
                st.markdown(f"**H-Bond Acceptors:** {props['hba_count']}")
                st.markdown(f"**Rotatable Bonds:** {props['rotatable_bonds']}")
                st.markdown(f"**Heavy Atoms:** {props['heavy_atoms']}")
                st.markdown(f"**Rings:** {props['rings']}")
                st.markdown(f"**Lipinski Violations:** {props['lipinski_violations']}")
    
    with col2:
        st.subheader("🧪 Molecular Structure")
        
        # 2D Structure
        st.markdown("**2D Structure:**")
        svg_2d = render_molecule_2d(molecule)
        if svg_2d:
            st.markdown(svg_2d, unsafe_allow_html=True)
        else:
            st.info("2D structure not available")
        
        # 3D Structure
        st.markdown("**3D Interactive Structure:**")
        render_molecule_3d(molecule)  # Function handles display internally

def main():
    st.title("🧬 FoTChemistry Discovery Dashboard")
    st.markdown("""
    **Quantum-guided molecular discovery system** • Real chemistry, no duplicates • Field of Truth methodology
    """)
    
    # Sidebar
    with st.sidebar:
        st.header("⚡ System Status")
        
        # Feature availability - focus on what works
        st.markdown("**🚀 Core Features:**")
        st.markdown("✅ Molecular Discovery Data")
        st.markdown("✅ Property Analysis") 
        st.markdown("✅ Statistical Dashboard")
        st.markdown("✅ Interactive Visualization")
        
        st.markdown("**🧬 Advanced Features:**")
        st.markdown("✅ 2D Molecular Structures (Pre-generated)")
        
        # 3D visualization status
        if HAS_STMOL and HAS_PY3DMOL:
            st.markdown("✅ 3D Visualization (Stmol + Py3Dmol)")
        elif HAS_STMOL:
            st.markdown("✅ 3D Visualization (Stmol)")
        elif HAS_PY3DMOL:
            st.markdown("✅ 3D Visualization (Py3Dmol)")
        else:
            st.markdown("📊 3D Visualization (Static)")
            
        st.markdown(f"{'✅' if HAS_QUANTUM else '❌'} Quantum Engine")
        st.markdown(f"{'✅' if HAS_AKG else '📁'} Database {'(Live)' if HAS_AKG else '(Static)'}")
        
        # Clear cache button
        if st.button("🗑️ Clear Cache"):
            st.cache_data.clear()
            st.success("✅ Cache cleared!")
    
    # Load discovery data
    discovery_data = load_discovery_data()
    
    if not discovery_data:
        st.error("❌ No discovery data available. Run discovery campaigns to generate molecules.")
        st.markdown("""
        **🚀 To start discovering molecules:**
        
        ```bash
        python3 continuous_chemistry_discovery.py --continuous
        ```
        """)
        return
    
    # Statistics overview
    stats = discovery_data.get('statistics', {})
    
    st.subheader("📈 Discovery Overview")
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("🧬 Total Molecules", stats.get('total_molecules', 0))
    with col2:
        st.metric("📊 Average Score", f"{stats.get('avg_score', 0):.3f}")
    with col3:
        st.metric("🏆 Max Score", f"{stats.get('max_score', 0):.3f}")
    with col4:
        st.metric("🔬 Active Claims", stats.get('active_claims', 0))
    
    # Main discovery list
    st.subheader("🧬 Discovered Molecules")
    
    # Use all discoveries, not just recent_molecules
    molecules = discovery_data.get('discoveries', discovery_data.get('recent_molecules', []))
    if not molecules:
        st.warning("⚠️ No molecules found in discovery data")
        return
    
    # Create a selectbox for molecule navigation
    molecule_options = [f"{i+1}. {mol.get('smiles', 'Unknown')} (Score: {mol.get('score', 0):.3f})" 
                       for i, mol in enumerate(molecules)]
    
    selected_idx = st.selectbox(
        "🔍 Select a molecule to view in detail:",
        range(len(molecule_options)),
        format_func=lambda x: molecule_options[x],
        index=0
    )
    
    # Display selected molecule
    if selected_idx is not None:
        selected_molecule = molecules[selected_idx]
        
        st.markdown("---")
        display_molecule_detail(selected_molecule)
    
    # Quick overview table
    st.markdown("---")
    st.subheader("📋 All Discoveries - Quick View")
    
    # Create summary dataframe
    summary_data = []
    for i, mol in enumerate(molecules):
        smiles = mol.get('smiles', '')
        summary_data.append({
            '#': i + 1,
            'SMILES': smiles if len(smiles) <= 25 else smiles[:22] + '...',
            'Score': f"{mol.get('score', 0):.3f}",
            'Drug Like': '✅' if (isinstance(mol.get('drug_likeness'), dict) and 
                                mol.get('drug_likeness', {}).get('passes_lipinski', False)) else '❌',
            'Safety': f"{mol.get('safety_score', 0):.3f}",
            'Coherence': f"{mol.get('quantum_coherence', 0):.6f}",
            'ID': mol.get('id', 'Unknown')[:8] + '...' if mol.get('id') else 'Unknown'
        })
    
    summary_df = pd.DataFrame(summary_data)
    st.dataframe(summary_df, use_container_width=True)
    
    # Footer info
    st.markdown("---")
    st.markdown("""
    **🎯 Discovery Quality:**
    - All molecules are unique (no duplicates)
    - Minimum complexity requirements enforced
    - Real chemical structures only
    - Quantum-guided property optimization
    
    **🔄 To refresh data:** Clear cache in sidebar and reload page
    """)

if __name__ == "__main__":
    main()
