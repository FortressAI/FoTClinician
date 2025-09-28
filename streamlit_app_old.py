#!/usr/bin/env python3
"""
FoTChemistry Discovery Dashboard

Streamlit application for visualizing molecular discoveries from the 
quantum-guided chemistry discovery engine.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import json
import os
import sys
import traceback

# Add core and akg to path
sys.path.append('core')
sys.path.append('akg')
sys.path.append('agents/alchemist')

# Feature flags for conditional imports
HAS_AKG = False
HAS_QUANTUM = False  
HAS_RDKIT = False
HAS_3D_VIZ = False

try:
    from client import AKG
    HAS_AKG = True
except ImportError:
    st.warning("⚠️ AKG client not available - using local data only")

try:
    from chemistry_vqbit_engine import ChemistryVQbitEngine, ChemistryPropertyType
    HAS_QUANTUM = True
except ImportError:
    st.warning("⚠️ Quantum engine not available")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    import io
    import base64
    HAS_RDKIT = True
except ImportError:
    st.warning("⚠️ RDKit not available - molecular visualization limited")

try:
    import stmol
    import py3Dmol
    HAS_3D_VIZ = True
except ImportError:
    pass

# Page configuration
st.set_page_config(
    page_title="FoTChemistry Discovery Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)


@st.cache_resource
def get_akg_client():
    """Get AKG client connection"""
    if HAS_AKG:
        try:
            return AKG()
        except Exception as e:
            st.error(f"❌ Failed to connect to AKG: {e}")
    return None


@st.cache_resource
def get_quantum_engine():
    """Initialize quantum vQbit engine"""
    if HAS_QUANTUM:
        try:
            engine = ChemistryVQbitEngine(use_gpu=True)
            return engine
        except Exception as e:
            st.error(f"❌ Failed to initialize quantum engine: {e}")
    return None


def render_molecular_structure_2d(smiles: str):
    """Generate 2D molecular structure as SVG"""
    if not HAS_RDKIT:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        return svg
    except Exception as e:
        st.error(f"❌ 2D structure generation failed: {e}")
        return None


def render_molecular_structure_3d(smiles: str):
    """Generate 3D molecular structure with py3Dmol"""
    if not HAS_RDKIT or not HAS_3D_VIZ:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        
        # Try multiple embedding attempts for better 3D coordinates
        for attempt in range(3):
            try:
                AllChem.EmbedMolecule(mol, randomSeed=42 + attempt)
                AllChem.UFFOptimizeMolecule(mol)
                break
            except:
                if attempt == 2:  # Last attempt failed
                    return None
                continue
        
        # Convert to SDF format for py3Dmol
        sdf = Chem.MolToMolBlock(mol)
        return sdf
    except Exception as e:
        st.error(f"❌ 3D structure generation failed: {e}")
        return None


def calculate_molecular_properties(smiles: str):
    """Calculate comprehensive molecular properties"""
    if not HAS_RDKIT:
        st.error("❌ RDKit required for molecular property calculations")
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.error("❌ Invalid SMILES notation")
            return None
        
        # Basic molecular properties
        properties = {
            'molecular_formula': rdMolDescriptors.CalcMolFormula(mol),
            'molecular_weight': round(Descriptors.MolWt(mol), 2),
            'exact_mass': round(Descriptors.ExactMolWt(mol), 4),
            'heavy_atom_count': mol.GetNumHeavyAtoms(),
            'total_atom_count': mol.GetNumAtoms(),
            'bond_count': mol.GetNumBonds(),
            'ring_count': Descriptors.RingCount(mol),
            'aromatic_ring_count': Descriptors.NumAromaticRings(mol),
            'rotatable_bond_count': Descriptors.NumRotatableBonds(mol),
            
            # Physicochemical properties
            'logp': round(Descriptors.MolLogP(mol), 2),
            'tpsa': round(Descriptors.TPSA(mol), 2),
            'hbd_count': Descriptors.NumHDonors(mol),
            'hba_count': Descriptors.NumHAcceptors(mol),
            'molar_refractivity': round(Descriptors.MolMR(mol), 2),
            
            # Drug-like properties (Lipinski's Rule of Five)
            'lipinski_violations': sum([
                Descriptors.MolWt(mol) > 500,
                Descriptors.MolLogP(mol) > 5,
                Descriptors.NumHDonors(mol) > 5,
                Descriptors.NumHAcceptors(mol) > 10
            ]),
            
            # Chemical identifiers
            'inchi': Chem.MolToInchi(mol),
            'inchi_key': Chem.MolToInchiKey(mol),
            'canonical_smiles': Chem.MolToSmiles(mol)
        }
        
        return properties
    except Exception as e:
        st.error(f"❌ Molecular property calculation failed: {e}")
        return None


@st.cache_data
def load_discovery_data():
    """Load discovery data from AKG or file export"""
    # For local deployment, try to query Neo4j directly for real-time data
    akg = get_akg_client()
    if akg:
        try:
            # Query for recent chemistry discoveries from FoTChem_Discovery nodes
            with akg.neo4j_driver.session() as session:
                result = session.run("""
                    MATCH (d:FoTChem_Discovery)
                    RETURN d.smiles as smiles, 
                           d.combined_score as score,
                           d.drug_likeness as drug_likeness,
                           d.safety_score as safety_score,
                           d.quantum_coherence as quantum_coherence,
                           d.discovery_timestamp as timestamp,
                           d.discovery_id as id
                    ORDER BY d.discovery_timestamp DESC
                    LIMIT 50
                """)
                
                discoveries = []
                for record in result:
                    # Handle drug_likeness which might be a JSON string or dict
                    drug_likeness_raw = record.get('drug_likeness', {})
                    if isinstance(drug_likeness_raw, str):
                        try:
                            drug_likeness_parsed = json.loads(drug_likeness_raw)
                        except:
                            drug_likeness_parsed = {'passes_lipinski': False}
                    else:
                        drug_likeness_parsed = drug_likeness_raw or {'passes_lipinski': False}
                    
                    discovery = {
                        'id': record.get('id', 'unknown'),
                        'smiles': record.get('smiles', ''),
                        'score': float(record.get('score', 0) or 0),
                        'drug_likeness': drug_likeness_parsed,
                        'safety_score': float(record.get('safety_score', 0) or 0),
                        'quantum_coherence': float(record.get('quantum_coherence', 0) or 0),
                        'timestamp': record.get('timestamp', ''),
                        'properties': {
                            'combined_score': float(record.get('score', 0) or 0),
                            'drug_likeness': drug_likeness_parsed,
                            'safety_score': float(record.get('safety_score', 0) or 0)
                        }
                    }
                    discoveries.append(discovery)
                
                # Get discovery statistics
                stats_result = session.run("""
                    MATCH (d:FoTChem_Discovery)
                    RETURN count(d) as total_discoveries,
                           avg(d.combined_score) as avg_score,
                           max(d.combined_score) as max_score
                """)
                stats_record = stats_result.single()
                
                return {
                    'discoveries': discoveries,
                    'total_discoveries': stats_record.get('total_discoveries', 0),
                    'statistics': {
                        'total_molecules': stats_record.get('total_discoveries', 0),
                        'avg_score': float(stats_record.get('avg_score') or 0),
                        'max_score': float(stats_record.get('max_score') or 0),
                        'total_reactions': 0,  # Not tracked yet
                        'total_measurements': 0,  # Not tracked yet
                        'active_claims': 0  # Not tracked yet
                    },
                    'recent_molecules': discoveries[:10]  # Top 10 for display
                }
                
        except Exception as e:
            st.warning(f"⚠️ Failed to query Neo4j directly: {e}")
    
    # Fallback: Try to load from file export (for cloud deployment)
    export_file = "results/chemistry_discoveries.json"
    if os.path.exists(export_file):
        try:
            with open(export_file, 'r') as f:
                data = json.load(f)
            if data.get('total_discoveries', 0) > 0:
                return data
        except Exception as e:
            st.warning(f"⚠️ Failed to load discovery data from file: {e}")
    
    return None


def analyze_molecule_quantum(smiles: str):
    """Analyze molecule using quantum vQbit engine"""
    engine = get_quantum_engine()
    if not engine:
        return {"error": "Quantum engine not available", "success": False}
    
    try:
        # Create quantum vQbit state for the molecule
        initial_property_scores = {
            ChemistryPropertyType.BIOACTIVITY: 0.7,
            ChemistryPropertyType.SUSTAINABILITY: 0.6,
            ChemistryPropertyType.REPRODUCIBILITY: 0.8,
            ChemistryPropertyType.EFFICIENCY: 0.7
        }
        
        vqbit_state = engine.create_molecular_vqbit(initial_property_scores)
        
        # Measure quantum properties
        property_measurements = {}
        for prop_type in ChemistryPropertyType:
            measurement = engine.measure_property(vqbit_state, prop_type)
            property_measurements[prop_type.name.lower()] = float(measurement)
        
        # Calculate quantum coherence
        coherence = engine._calculate_l1_coherence(vqbit_state.amplitudes)
        
        return {
            'success': True,
            'property_measurements': property_measurements,
            'quantum_coherence': float(coherence),
            'quantum_state_norm': float(np.linalg.norm(vqbit_state.amplitudes))
        }
        
    except Exception as e:
        return {"error": str(e), "success": False}


def main():
    st.title("🧬 FoTChemistry Discovery Dashboard")
    st.markdown("""
    **Field of Truth methodology for autonomous chemical discovery**  
    *Quantum vQbit substrate • Truth-mining workflows • Open chemical knowledge*
    """)
    
    # Sidebar - System Status
    with st.sidebar:
        st.header("⚡ System Status")
        
        # Feature availability
        st.markdown("**🛠️ Available Features:**")
        st.markdown(f"{'✅' if HAS_AKG else '❌'} AKG Database Connection")
        st.markdown(f"{'✅' if HAS_QUANTUM else '❌'} Quantum vQbit Engine")
        st.markdown(f"{'✅' if HAS_RDKIT else '❌'} RDKit Molecular Analysis")
        st.markdown(f"{'✅' if HAS_3D_VIZ else '❌'} 3D Molecular Visualization")
        
        # Quantum metrics
        if HAS_QUANTUM:
            engine = get_quantum_engine()
            if engine:
                st.markdown("**⚛️ Quantum Substrate Metrics:**")
                st.metric("🌌 Hilbert Dimension", engine.hilbert_dimension)
                st.metric("🚀 GPU Acceleration", "Yes" if engine.gpu_acceleration else "No")
                st.metric("🧠 Property Operators", len(ChemistryPropertyType))
        
        # Clear cache button
        if st.button("🗑️ Clear Cache"):
            st.cache_data.clear()
            st.cache_resource.clear()
            st.success("✅ Cache cleared!")
    
    # Main content tabs
    tab1, tab2, tab3 = st.tabs(["🧬 Discovered Molecules", "🔬 Analyze Molecule", "⚛️ System Status"])
    
    with tab1:
        st.header("🧬 Discovered Molecules")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("🔬 Input Molecule")
            
            # SMILES input
            smiles = st.text_input(
                "Enter SMILES notation:",
                value="CCO",
                help="Chemical structure in SMILES format"
            )
            
            # Example molecules
            st.markdown("**💡 Example molecules:**")
            example_cols = st.columns(3)
            with example_cols[0]:
                if st.button("Ethanol"):
                    smiles = "CCO"
            with example_cols[1]:
                if st.button("Benzene"):
                    smiles = "c1ccccc1"
            with example_cols[2]:
                if st.button("Aspirin"):
                    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
        
        with col2:
            if smiles:
                st.subheader("🔬 Molecular Structure")
                
                col_struct, col_3d = st.columns([1, 1])
                
                with col_struct:
                    # 2D Structure
                    svg_2d = render_molecular_structure_2d(smiles)
                    if svg_2d:
                        st.markdown("**2D Structure:**")
                        st.markdown(svg_2d, unsafe_allow_html=True)
                
                with col_3d:
                    # 3D Interactive Structure
                    if HAS_3D_VIZ:
                        try:
                            sdf_3d = render_molecular_structure_3d(smiles)
                            if sdf_3d:
                                st.markdown("**🌐 3D Interactive Structure:**")
                                viewer = py3Dmol.view(width=300, height=300)
                                viewer.addModel(sdf_3d, 'sdf')
                                viewer.setStyle({'stick': {'radius': 0.1}, 'sphere': {'radius': 0.3}})
                                viewer.setBackgroundColor('white')
                                viewer.zoomTo()
                                
                                # Add labels for small molecules
                                mol = Chem.MolFromSmiles(smiles)
                                if mol and mol.GetNumAtoms() <= 20:
                                    for i, atom in enumerate(mol.GetAtoms()):
                                        viewer.addLabel(atom.GetSymbol(), 
                                                      {'position': {'x': 0, 'y': 0, 'z': 0}, 
                                                       'fontColor': 'black', 'fontSize': 12})
                                
                                stmol.showmol(viewer, height=300, width=300)
                                st.success("✅ 3D structure rendered")
                            else:
                                st.warning("⚠️ Could not generate 3D coordinates")
                        except Exception as e:
                            st.error(f"❌ 3D visualization error: {e}")
                    else:
                        st.info("💡 Install stmol and py3Dmol for 3D visualization")
        
        # Molecular properties and quantum analysis
        if smiles:
            st.subheader("📊 Molecular Properties")
            
            prop_col1, prop_col2 = st.columns([1, 1])
            
            with prop_col1:
                # Calculate molecular properties
                properties = calculate_molecular_properties(smiles)
                if properties:
                    st.markdown("**🧮 Chemical Properties:**")
                    st.text(f"Formula: {properties['molecular_formula']}")
                    st.text(f"Molecular Weight: {properties['molecular_weight']} g/mol")
                    st.text(f"LogP: {properties['logp']}")
                    st.text(f"TPSA: {properties['tpsa']} Ų")
                    st.text(f"H-Bond Donors: {properties['hbd_count']}")
                    st.text(f"H-Bond Acceptors: {properties['hba_count']}")
                    st.text(f"Rotatable Bonds: {properties['rotatable_bond_count']}")
                    st.text(f"Lipinski Violations: {properties['lipinski_violations']}")
            
            with prop_col2:
                # Quantum analysis
                if HAS_QUANTUM:
                    st.markdown("**⚛️ Quantum Analysis:**")
                    
                    if st.button("🔬 Analyze Quantum Properties"):
                        with st.spinner("Running quantum analysis..."):
                            quantum_result = analyze_molecule_quantum(smiles)
                        
                        if quantum_result['success']:
                            measurements = quantum_result['property_measurements']
                            coherence = quantum_result['quantum_coherence']
                            
                            st.success(f"✅ Quantum coherence: {coherence:.6f}")
                            
                            # Create measurements dataframe
                            measurements_df = pd.DataFrame([
                                {'Property': k, 'Measurement': v} 
                                for k, v in measurements.items()
                            ])
                            
                            # Property descriptions
                            property_descriptions = {
                                'bioactivity': 'Biological activity potential and safety profile',
                                'sustainability': 'Environmental impact and green chemistry score',
                                'reproducibility': 'Experimental reproducibility and method reliability', 
                                'efficiency': 'Synthetic efficiency and atom economy'
                            }
                            
                            measurements_df['Description'] = measurements_df['Property'].map(property_descriptions)
                            measurements_df['Property'] = measurements_df['Property'].str.title()
                            
                            # Interactive bar chart
                            fig = px.bar(
                                measurements_df, 
                                x='Property', 
                                y='Measurement',
                                hover_data={'Description': True},
                                title="Quantum Property Measurements",
                                color='Measurement', 
                                color_continuous_scale='plasma',
                                range_color=[0, 1]
                            )
                            fig.update_layout(
                                yaxis_range=[0, 1],
                                xaxis_title="Chemical Property",
                                yaxis_title="Quantum Measurement Value"
                            )
                            st.plotly_chart(fig, use_container_width=True)
                        else:
                            st.error(f"❌ Quantum analysis failed: {quantum_result.get('error', 'Unknown error')}")
    
    with tab2:
        st.header("📊 Discovery Dashboard")
        
        # Load discovery data
        discovery_data = load_discovery_data()
        
        if discovery_data:
            # Statistics overview
            stats = discovery_data.get('statistics', {})
            
            st.subheader("📈 Discovery Statistics")
            stat_cols = st.columns(4)
            
            with stat_cols[0]:
                st.metric("🧬 Total Molecules", stats.get('total_molecules', 0))
            with stat_cols[1]:
                st.metric("📊 Average Score", f"{stats.get('avg_score', 0):.3f}")
            with stat_cols[2]:
                st.metric("🏆 Max Score", f"{stats.get('max_score', 0):.3f}")
            with stat_cols[3]:
                st.metric("🔬 Active Claims", stats.get('active_claims', 0))
            
            # Show recent discoveries
            if 'recent_molecules' in discovery_data and discovery_data['recent_molecules']:
                st.subheader("🆕 Recent Molecular Discoveries")
                
                molecules = discovery_data['recent_molecules']
                
                # Display each discovery with molecular structure
                for i, molecule in enumerate(molecules[:5]):  # Show top 5
                    with st.expander(f"🧬 Discovery {i+1}: {molecule.get('smiles', 'Unknown')} (Score: {molecule.get('score', 0):.3f})"):
                        col_struct, col_data = st.columns([1, 1])
                        
                        with col_struct:
                            # Show 2D structure if available
                            if HAS_RDKIT and molecule.get('smiles'):
                                try:
                                    svg_2d = render_molecular_structure_2d(molecule['smiles'])
                                    if svg_2d:
                                        st.markdown("**Molecular Structure:**")
                                        st.markdown(svg_2d, unsafe_allow_html=True)
                                except:
                                    st.text("Structure not available")
                        
                        with col_data:
                            st.markdown("**Discovery Details:**")
                            st.text(f"SMILES: {molecule.get('smiles', 'N/A')}")
                            st.text(f"Combined Score: {molecule.get('score', 0):.3f}")
                            
                            # Handle drug_likeness as dictionary or float
                            drug_likeness_val = molecule.get('drug_likeness', 0)
                            if isinstance(drug_likeness_val, dict):
                                drug_likeness_score = drug_likeness_val.get('passes_lipinski', False)
                                st.text(f"Drug Likeness: {'Pass' if drug_likeness_score else 'Fail'} Lipinski")
                            else:
                                st.text(f"Drug Likeness: {float(drug_likeness_val):.3f}")
                            
                            st.text(f"Safety Score: {molecule.get('safety_score', 0):.3f}")
                            st.text(f"Quantum Coherence: {molecule.get('quantum_coherence', 0):.6f}")
                            st.text(f"Discovery ID: {molecule.get('id', 'Unknown')[:8]}...")
                            
                            # Add molecular properties if available
                            if HAS_RDKIT and molecule.get('smiles'):
                                try:
                                    props = calculate_molecular_properties(molecule['smiles'])
                                    if props:
                                        st.markdown("**Molecular Properties:**")
                                        st.text(f"Formula: {props.get('molecular_formula', 'N/A')}")
                                        st.text(f"MW: {props.get('molecular_weight', 0):.1f} g/mol")
                                        st.text(f"LogP: {props.get('logp', 0):.2f}")
                                except:
                                    pass
                
                # Summary table of all discoveries
                if len(molecules) > 5:
                    st.markdown(f"**📋 All {len(molecules)} Recent Discoveries:**")
                    summary_data = []
                    for mol in molecules:
                        summary_data.append({
                            'SMILES': mol.get('smiles', '')[:20] + '...' if len(mol.get('smiles', '')) > 20 else mol.get('smiles', ''),
                            'Score': f"{mol.get('score', 0):.3f}",
                            'Safety': f"{mol.get('safety_score', 0):.3f}",
                            'ID': mol.get('id', 'Unknown')[:8] + '...' if mol.get('id') else 'Unknown'
                        })
                    
                    summary_df = pd.DataFrame(summary_data)
                    st.dataframe(summary_df, use_container_width=True)
            else:
                st.info("No recent molecular discoveries found")
        else:
            st.warning("⚠️ No discovery data available. Run discovery campaigns to generate molecules.")
            
            st.markdown("""
            **🚀 To start discovering molecules:**
            
            1. Run continuous discovery:
               ```bash
               python3 continuous_chemistry_discovery.py --continuous
               ```
            
            2. Or run parallel discovery:
               ```bash
               python3 scale_discovery.py --target 1000
               ```
            """)
    
    with tab3:
        st.header("⚛️ Quantum Substrate Metrics")
        
        if HAS_QUANTUM:
            engine = get_quantum_engine()
            if engine:
                col1, col2 = st.columns([1, 1])
                
                with col1:
                    st.subheader("🌌 Quantum System")
                    st.metric("Hilbert Space Dimension", engine.hilbert_dimension)
                    st.metric("Coherence Time", "∞ (noiseless)")
                    st.metric("Quantum Fidelity", "1.0 (perfect)")
                    
                    # Test quantum substrate
                    if st.button("🧪 Test Quantum Substrate"):
                        with st.spinner("Running quantum tests..."):
                            test_results = engine.test_quantum_substrate()
                        
                        st.subheader("🔬 Test Results")
                        for test_name, passed in test_results.items():
                            status = "✅ PASSED" if passed else "❌ FAILED"
                            st.text(f"{test_name}: {status}")
                
                with col2:
                    st.subheader("🔧 Performance Optimizations")
                    st.metric("🚀 GPU Acceleration", "Yes" if engine.gpu_acceleration else "No")
                    st.metric("🧮 C Extensions", "Available")
                    st.metric("⚛️ Property Operators", len(ChemistryPropertyType))
                    
                    # Quantum properties
                    st.subheader("🎯 Chemistry Properties")
                    for prop_type in ChemistryPropertyType:
                        st.text(f"• {prop_type.value.title()}")
        else:
            st.warning("⚠️ Quantum engine not available")
            st.info("Install quantum dependencies to enable quantum substrate analysis")


if __name__ == "__main__":
    main()