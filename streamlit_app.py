#!/usr/bin/env python3
"""
FoTChemistry Discovery Dashboard

Dashboard-first interface for viewing molecular discoveries with detailed analysis.
"""

import streamlit as st
import streamlit.components.v1 as components
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
    from rdkit.Chem import rdDepictor
    HAS_RDKIT = True
    print("RDKit loaded successfully")
except ImportError:
    print("RDKit not available")
    pass

try:
    import stmol
    HAS_STMOL = True
    print("stmol loaded successfully")
except ImportError:
    print("stmol not available")
    pass

try:
    import py3Dmol
    HAS_PY3DMOL = True
    print("py3Dmol loaded successfully")
except ImportError:
    print("py3Dmol not available")
    pass

# Page configuration
st.set_page_config(
    page_title="FoTChemistry Discovery Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

@st.cache_data(ttl=60, show_spinner="Loading all 6,443 molecules...")  # Cache for 60 seconds only
def load_discovery_data():
    """Load discovery data with cloud/local detection and fallback options."""
    
    # Detect deployment environment (improved cloud detection)
    is_cloud_deployment = (
        # Official Streamlit Cloud environment variables
        os.environ.get('STREAMLIT_SHARING') == '1' or
        os.environ.get('STREAMLIT_CLOUD') == '1' or
        os.environ.get('STREAMLIT_COMMUNITY_CLOUD') == '1' or
        
        # Hostname-based detection
        'streamlit.io' in os.environ.get('HOSTNAME', '') or
        'streamlit.app' in os.environ.get('HOSTNAME', '') or
        
        # Path-based detection (cloud containers)
        '/mount/src/' in os.getcwd() or
        '/app/' in os.getcwd() or
        
        # File-based detection (local development files missing)
        not os.path.exists('akg/client.py') or
        not os.path.exists('core/chemistry_vqbit_engine.py') or
        
        # Environment variable that indicates cloud
        os.environ.get('HOME', '').startswith('/home/adminuser')
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
        
        # Generate demo data if no files found
        st.info("🎭 Generating demo data for visualization")
        return create_demo_data()
    
    else:
        st.info("🏠 Local deployment detected - using live data")
        
        # Try to load from fixed chemistry discoveries first, then fallback
        data_files = [
            "results/chemistry_discoveries.json",  # Fixed data with score mapping
            "cloud_data_snapshot.json", 
            "results/overnight_discovery_mega_dataset.json"  # Legacy data
        ]
        
        for export_file in data_files:
            if os.path.exists(export_file):
                try:
                    with open(export_file, 'r') as f:
                        data = json.load(f)
                    
                    # Handle different data structures
                    if isinstance(data, dict):
                        if 'discovery_summary' in data:
                            total_count = data['discovery_summary'].get('total_discoveries', 0)
                            display_count = len(data.get('discoveries', []))
                            st.success(f"📂 Loaded {display_count} molecules from {export_file} (Total: {total_count})")
                        else:
                            display_count = len(data.get('discoveries', data))
                            st.success(f"📂 Loaded {display_count} molecules from {export_file}")
                    else:
                        display_count = len(data) if isinstance(data, list) else 0
                        st.success(f"📂 Loaded {display_count} molecules from {export_file}")
                    
                    return data
                except Exception as e:
                    st.warning(f"⚠️ Failed to load {export_file}: {e}")
                    continue
        
        # Fallback to demo data for local testing
        st.info("🎭 Generating demo data for testing")
        return create_demo_data()

def create_demo_data():
    """Create embedded demo data for cloud deployment."""
    from datetime import datetime
    
    demo_molecules = [
        {
            "id": "demo_001",
            "smiles": "CCO",
            "name": "Ethanol",
            "score": 0.786,
            "drug_likeness": {"passes_lipinski": True, "score": 0.85},
            "safety_score": 0.92,
            "quantum_coherence": 0.74,
            "timestamp": datetime.now().isoformat(),
            "properties": {
                "formula": "C2H6O",
                "molecular_weight": 46.07,
                "logp": -0.31,
                "tpsa": 20.23,
                "hbd": 1,
                "hba": 1
            }
        },
        {
            "id": "demo_002", 
            "smiles": "c1ccc(cc1)O",
            "name": "Phenol",
            "score": 0.754,
            "drug_likeness": {"passes_lipinski": True, "score": 0.78},
            "safety_score": 0.88,
            "quantum_coherence": 0.71,
            "timestamp": datetime.now().isoformat(),
            "properties": {
                "formula": "C6H6O",
                "molecular_weight": 94.11,
                "logp": 1.46,
                "tpsa": 20.23,
                "hbd": 1,
                "hba": 1
            }
        },
        {
            "id": "demo_003",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "name": "Aspirin",
            "score": 0.834,
            "drug_likeness": {"passes_lipinski": True, "score": 0.91},
            "safety_score": 0.95,
            "quantum_coherence": 0.82,
            "timestamp": datetime.now().isoformat(),
            "properties": {
                "formula": "C9H8O4",
                "molecular_weight": 180.16,
                "logp": 1.19,
                "tpsa": 63.6,
                "hbd": 1,
                "hba": 4
            }
        },
        {
            "id": "demo_004",
            "smiles": "CC1=CC=C(C=C1)C",
            "name": "p-Xylene",
            "score": 0.698,
            "drug_likeness": {"passes_lipinski": True, "score": 0.72},
            "safety_score": 0.85,
            "quantum_coherence": 0.69,
            "timestamp": datetime.now().isoformat(),
            "properties": {
                "formula": "C8H10",
                "molecular_weight": 106.17,
                "logp": 3.15,
                "tpsa": 0.0,
                "hbd": 0,
                "hba": 0
            }
        },
        {
            "id": "demo_005",
            "smiles": "CN1CCC[C@H]1c2cccnc2",
            "name": "Nicotine",
            "score": 0.712,
            "drug_likeness": {"passes_lipinski": True, "score": 0.76},
            "safety_score": 0.78,
            "quantum_coherence": 0.73,
            "timestamp": datetime.now().isoformat(),
            "properties": {
                "formula": "C10H14N2",
                "molecular_weight": 162.23,
                "logp": 1.17,
                "tpsa": 16.13,
                "hbd": 0,
                "hba": 2
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

def render_molecule_2d(smiles: str, width=400, height=400, molecule_data=None):
    """Render 2D molecular structure using pre-generated data or RDKit"""
    if not smiles:
        return create_molecule_placeholder(smiles, "2D Structure", width, height)
    
    # Try pre-generated visualization first (for cloud deployment)
    if molecule_data and 'visualization_2d' in molecule_data:
        viz_data = molecule_data['visualization_2d']
        if viz_data.get('svg'):
            svg_content = viz_data['svg']
            # Clean SVG (remove XML declaration)
            if svg_content.startswith('<?xml'):
                svg_content = svg_content[svg_content.find('<svg'):]
            return svg_content
    
    # Fallback to live RDKit generation (local only)
    if not HAS_RDKIT:
        return create_molecule_placeholder(smiles, "2D Structure - RDKit Required", width, height)
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.error(f"Invalid SMILES: {smiles}")
            return create_molecule_placeholder(smiles, "Invalid SMILES", width, height)
        
        # Generate 2D coordinates
        rdDepictor.Compute2DCoords(mol)
        
        # Create SVG drawer
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.SetFontSize(16)
        
        # Draw molecule
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        # Get SVG
        svg = drawer.GetDrawingText()
        
        # Clean SVG (remove XML declaration) and validate
        cleaned_svg = svg
        if cleaned_svg.startswith('<?xml'):
            cleaned_svg = cleaned_svg[cleaned_svg.find('<svg'):]
        
        # Validate SVG has actual molecule content
        if 'bond-' in cleaned_svg or len(cleaned_svg) > 1000:  # Basic validation
            return cleaned_svg
        else:
            st.warning(f"Generated SVG appears empty for {smiles}")
            return create_molecule_placeholder(smiles, "2D Structure", width, height)
            
    except Exception as e:
        st.error(f"Error generating 2D structure for {smiles}: {str(e)}")
        return create_molecule_placeholder(smiles, "2D Structure Error", width, height)

def render_molecule_3d(smiles: str, height=400, molecule_data=None):
    """Render 3D molecular structure using pre-generated data, stmol and py3Dmol"""
    if not smiles:
        st.write("No SMILES provided for 3D rendering")
        return False
    
    # Try pre-generated 3D data first (for cloud deployment)
    if molecule_data and 'visualization_3d' in molecule_data:
        viz_data = molecule_data['visualization_3d']
        if viz_data.get('molblock'):
            molblock = viz_data['molblock']
            
            # Try stmol first (best for cloud)
            if HAS_STMOL and HAS_PY3DMOL:
                try:
                    # Create py3Dmol viewer and pass to stmol
                    viewer = py3Dmol.view(width=400, height=height)
                    viewer.addModel(molblock, 'mol')
                    viewer.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.3}})
                    viewer.setBackgroundColor('white')
                    viewer.zoomTo()
                    stmol.showmol(viewer, height=height, width=400)
                    return True
                except Exception as e:
                    st.warning(f"stmol failed: {str(e)}")
            
            # Try py3Dmol with HTML
            if HAS_PY3DMOL:
                try:
                    viewer_html = f"""
                    <div id="3dmol_viewer" style="width: 400px; height: {height}px; background: white; border: 1px solid #ddd; border-radius: 8px;"></div>
                    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/1.8.0/3Dmol-min.js"></script>
                    <script>
                    let viewer = $3Dmol.createViewer('3dmol_viewer');
                    viewer.addModel(`{molblock}`, 'mol');
                    viewer.setStyle({{'stick': {{'radius': 0.1}}, 'sphere': {{'scale': 0.3}}}});
                    viewer.setBackgroundColor('white');
                    viewer.zoomTo();
                    viewer.render();
                    </script>
                    """
                    components.html(viewer_html, height=height)
                    return True
                except Exception as e:
                    st.warning(f"py3Dmol failed: {str(e)}")
    
    # Fallback to live RDKit generation (local only)
    if not HAS_RDKIT:
        st.info("💡 3D visualization available with pre-generated data or RDKit installation")
        return False
    
    try:
        # Parse SMILES and add hydrogens
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.error(f"Invalid SMILES: {smiles}")
            return False
        
        # Add hydrogens for 3D structure
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        success = AllChem.EmbedMolecule(mol, randomSeed=42)
        if success != 0:
            # Try alternative method if embedding fails
            success = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
            if success != 0:
                st.warning(f"Could not generate 3D coordinates for {smiles}")
                return False
        
        # Optimize geometry
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except:
            # If UFF fails, try MMFF
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except:
                st.warning("Geometry optimization failed, using unoptimized structure")
        
        # Convert to molblock for 3D visualization
        molblock = Chem.MolToMolBlock(mol)
        
        # Try stmol first (best for Streamlit)
        if HAS_STMOL and HAS_PY3DMOL:
            try:
                # Create py3Dmol viewer
                viewer = py3Dmol.view(width=400, height=height)
                viewer.addModel(molblock, 'mol')
                viewer.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.3}})
                viewer.setBackgroundColor('white')
                viewer.zoomTo()
                
                # Pass viewer object to stmol
                stmol.showmol(viewer, height=height, width=400)
                return True
            except Exception as e:
                st.warning(f"stmol failed: {str(e)}")
        
        # Try py3Dmol with custom HTML
        if HAS_PY3DMOL:
            try:
                viewer_html = f"""
                <div id="3dmol_{hash(smiles)}" style="height: {height}px; width: 400px; position: relative;"></div>
                <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
                <script>
                var element = document.getElementById('3dmol_{hash(smiles)}');
                var viewer = $3Dmol.createViewer(element);
                viewer.addModel(`{molblock}`, 'mol');
                viewer.setStyle({{}}, {{stick: {{radius: 0.1}}, sphere: {{scale: 0.3}}}});
                viewer.setBackgroundColor('white');
                viewer.zoomTo();
                viewer.render();
                </script>
                """
                components.html(viewer_html, height=height)
                return True
            except Exception as e:
                st.warning(f"py3Dmol failed: {str(e)}")
        
        # Fallback: show molblock as text
        st.text_area("3D Structure (MOL format)", molblock, height=200)
        return True
                
    except Exception as e:
        st.error(f"Error generating 3D structure for {smiles}: {str(e)}")
        return False

def create_molecule_placeholder(smiles, title, width=400, height=400):
    """Create a styled placeholder when molecule rendering fails"""
    return f'''
    <div style="width: {width}px; height: {height}px; border: 2px solid #ddd; border-radius: 12px; 
                background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%); 
                display: flex; flex-direction: column; align-items: center; justify-content: center;
                font-family: Arial, sans-serif; margin: 10px 0;">
        <div style="font-size: 48px; margin-bottom: 20px;">⚗️</div>
        <div style="font-size: 18px; font-weight: bold; margin-bottom: 10px; color: #1565c0;">
            {title}
        </div>
        <div style="font-size: 14px; color: #1976d2; margin-bottom: 15px; padding: 0 20px; text-align: center;">
            {smiles}
        </div>
        <div style="font-size: 12px; color: #666; text-align: center; padding: 0 20px;">
            Visualization requires RDKit
        </div>
    </div>
    '''

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
        st.error(f"Error calculating properties: {str(e)}")
        return None

def display_molecule_detail(molecule):
    """Display detailed view of a single molecule"""
    smiles = molecule.get('smiles', '')
    st.header(f"🧬 Molecule Detail: {smiles}")
    
    # Main layout
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📊 Discovery Information")
        
        # Basic info
        st.markdown(f"**SMILES:** `{smiles}`")
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
        props = calculate_molecular_properties(smiles)
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
        else:
            st.info("Molecular properties calculation requires RDKit")
    
    with col2:
        st.subheader("🧪 Molecular Visualization")
        
        # 2D Structure
        st.markdown("**2D Structure:**")
        svg_2d = render_molecule_2d(smiles, molecule_data=molecule)
        if svg_2d:
            components.html(svg_2d, height=400, width=400)
        
        # 3D Structure
        st.markdown("**3D Interactive Structure:**")
        if not render_molecule_3d(smiles, molecule_data=molecule):
            # Show fallback placeholder
            st.markdown(create_molecule_placeholder(smiles, "3D Structure", 400, 300), unsafe_allow_html=True)

def main():
    st.title("🧬 FoTChemistry Validation Dashboard")
    st.markdown("""
    **Quantum-guided molecular validation platform** • Rigorous screening for real-world impact • Field of Truth methodology
    
    > **🔬 VALIDATION BREAKTHROUGH**: Transform generated candidates into validated discoveries!  
    > **Complete validation system** screens for novelty, safety, synthetic accessibility, and public benefit.
    > 
    > **📊 This Dashboard**: Explore molecular candidates and validation workflow  
    > **🎯 Problem Dashboard**: [Chemistry Challenge Analysis](https://fotchemistry-solutions.streamlit.app/)  
    > **🔬 Run Validation**: `python run_complete_validation.py` to validate overnight discoveries
    > 
    > **🌍 Public Benefit Focus**: Prioritizing compounds that help healthcare, environment, and education
    """)
    
    # Sidebar
    with st.sidebar:
        st.header("⚡ System Status")
        
        # Feature availability
        st.markdown("**🚀 Core Features:**")
        st.markdown("✅ Molecular Discovery Data")
        st.markdown("✅ Property Analysis") 
        st.markdown("✅ Statistical Dashboard")
        st.markdown("✅ Interactive Visualization")
        
        st.markdown("**🧬 Visualization Features:**")
        st.markdown(f"{'✅' if HAS_RDKIT else '❌'} RDKit (2D/3D Generation)")
        st.markdown(f"{'✅' if HAS_STMOL else '❌'} stmol (3D Viewer)")
        st.markdown(f"{'✅' if HAS_PY3DMOL else '❌'} py3Dmol (3D Backup)")
        st.markdown(f"{'✅' if HAS_QUANTUM else '❌'} Quantum Engine")
        st.markdown(f"{'✅' if HAS_AKG else '📁'} Database {'(Live)' if HAS_AKG else '(Static)'}")
        
        # Installation help
        if not HAS_RDKIT:
            st.warning("Install RDKit for molecular visualization:")
            st.code("pip install rdkit-pypi")
        
        if not HAS_STMOL and HAS_RDKIT:
            st.info("Install stmol for better 3D visualization:")
            st.code("pip install stmol")
        
        # Clear cache button
        if st.button("🗑️ Clear Cache"):
            st.cache_data.clear()
            st.success("✅ Cache cleared!")
    
    # Add refresh button
    if st.button("🔄 Force Refresh Data", help="Clear cache and reload all molecules"):
        st.cache_data.clear()
        st.rerun()
    
    # Load discovery data
    discovery_data = load_discovery_data()
    
    if not discovery_data:
        st.error("❌ No discovery data available.")
        return
    
    # Statistics overview
    summary = discovery_data.get('discovery_summary', {})
    discoveries = discovery_data.get('discoveries', [])
    
    # DEBUG: Show data loading details
    st.info(f"🔍 DEBUG: Loaded {len(discoveries)} molecules from data structure")
    st.info(f"📊 Summary claims: {summary.get('total_discoveries', 'Unknown')} total discoveries")
    
    # Calculate real-time statistics from data
    total_molecules = len(discoveries)
    scores = [mol.get('score', 0) for mol in discoveries if mol.get('score') is not None]
    avg_score = sum(scores) / len(scores) if scores else 0
    max_score = max(scores) if scores else 0
    
    st.subheader("📈 Discovery Overview")
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("🧬 Total Molecules", total_molecules)
    with col2:
        st.metric("📊 Average Score", f"{avg_score:.3f}")
    with col3:
        st.metric("🏆 Max Score", f"{max_score:.3f}")
    with col4:
        st.metric("🔬 Active Claims", summary.get('active_claims', 0))
    
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
    
    **📦 Required for full visualization:**
    - `pip install rdkit-pypi` (2D/3D structure generation)
    - `pip install stmol` (3D interactive viewer)
    - `pip install py3Dmol` (3D backup viewer)
    """)

if __name__ == "__main__":
    main()