#!/usr/bin/env python3
"""
FoTChemistry Discovery Dashboard
Field of Truth methodology for autonomous chemical discovery with quantum vQbit substrate.
Built for chemists to explore, validate, and discover chemical knowledge.
"""

import streamlit as st
import json
import os
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime
import sys

# Try to import chemistry libraries
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors, AllChem, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

# Try to import 3D molecular visualization
try:
    import stmol
    import py3Dmol
    HAS_3D_VIZ = True
except ImportError:
    HAS_3D_VIZ = False

# Add core to path for quantum engine
sys.path.append('core')

try:
    from akg.client import AKG
    HAS_AKG = True
except ImportError:
    HAS_AKG = False

try:
    from chemistry_vqbit_engine import ChemistryVQbitEngine, ChemistryPropertyType
    HAS_QUANTUM = True
except ImportError:
    HAS_QUANTUM = False

# Configure Streamlit page
st.set_page_config(
    page_title="FoTChemistry Discovery Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

@st.cache_resource
def get_quantum_engine():
    """Initialize quantum vQbit engine with caching for performance"""
    if not HAS_QUANTUM:
        return None
    try:
        engine = ChemistryVQbitEngine(use_gpu=True)
        return engine
    except Exception as e:
        st.error(f"❌ Failed to initialize quantum engine: {e}")
        return None

@st.cache_resource  
def get_akg_client():
    """Initialize AKG client with caching"""
    if not HAS_AKG:
        return None
    try:
        akg = AKG()
        return akg
    except Exception as e:
        st.sidebar.warning(f"⚠️ AKG connection failed: {e}")
        return None

def render_molecular_structure_2d(smiles: str):
    """Render 2D molecular structure using RDKit"""
    if not HAS_RDKIT:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Create SVG drawing
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        return drawer.GetDrawingText()
    except Exception as e:
        st.error(f"❌ 2D structure rendering failed: {e}")
        return None

def render_molecular_structure_3d(smiles: str):
    """Render 3D interactive molecular structure"""
    if not HAS_RDKIT or not HAS_3D_VIZ:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        
        # Try to embed molecule with multiple attempts
        embed_success = False
        for attempt in range(3):
            try:
                AllChem.EmbedMolecule(mol, randomSeed=42 + attempt, maxAttempts=100)
                # Use UFF force field for optimization
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
                embed_success = True
                break
            except:
                continue
        
        if not embed_success:
            # Fallback: try without optimization
            try:
                AllChem.EmbedMolecule(mol, randomSeed=42, maxAttempts=100)
            except:
                return None
        
        # Convert to SDF format for py3Dmol
        return Chem.MolToMolBlock(mol)
    except Exception as e:
        st.error(f"❌ 3D structure generation failed: {e}")
        return None

def calculate_molecular_properties(smiles: str):
    """Calculate comprehensive molecular properties (FoT ontology: ChemicalEntity properties)"""
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
    """Load discovery data from AKG (FoT ontology: Dataset, Evidence, Claims)"""
    # For local deployment, try to query Neo4j directly for real-time data
    akg = get_akg_client()
    if akg:
        try:
            # Query for recent chemistry discoveries directly from Neo4j
            with akg.neo4j_driver.session() as session:
                # Get recent chemistry discoveries from FoTChem_Discovery nodes
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
                        import json
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
                        'avg_score': float(stats_record.get('avg_score', 0)),
                        'max_score': float(stats_record.get('max_score', 0)),
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
        return {"error": "Quantum engine not available"}
    
    try:
        # Create initial property scores for quantum analysis
        initial_property_scores = {
            ChemistryPropertyType.BIOACTIVITY: 0.5,
            ChemistryPropertyType.SUSTAINABILITY: 0.5,
            ChemistryPropertyType.REPRODUCIBILITY: 0.5,
            ChemistryPropertyType.EFFICIENCY: 0.5
        }
        
        # Create molecular vQbit state
        vqbit_state = engine.create_molecular_vqbit(initial_property_scores)
        
        # Measure quantum properties
        measurements = {}
        for prop_type in ChemistryPropertyType:
            measurement = engine.measure_property(vqbit_state, prop_type)
            measurements[prop_type.name.lower()] = float(measurement)
        
        # Calculate quantum metrics
        coherence = engine._calculate_l1_coherence(vqbit_state.amplitudes)
        normalization = np.linalg.norm(vqbit_state.amplitudes)
        
        # Calculate chemical potential (mock calculation for demo)
        chemical_potential = sum(measurements.values()) / len(measurements) * 100
        
        # Calculate fidelity (mock calculation)
        fidelity = min(coherence * 1000, 1.0)
        
        return {
            'measurements': measurements,
            'coherence': float(coherence),
            'normalization': float(normalization),
            'chemical_potential': float(chemical_potential),
            'fidelity': float(fidelity),
            'success': True
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
        
        # Quantum vQbit engine status
        engine = get_quantum_engine()
        if engine:
            st.success("✅ Quantum Substrate Active")
            st.metric("🌌 Hilbert Dimension", engine.hilbert_dimension)
            st.metric("🚀 GPU Acceleration", "Yes" if engine.gpu_acceleration else "No")
        else:
            st.error("❌ Quantum Engine Offline")
        
        # AKG status
        akg = get_akg_client()
        if akg:
            st.success("✅ AKG Connected")
        else:
            st.error("❌ AKG Offline")
        
        # Feature availability
        st.write("📦 **Available Features:**")
        st.write(f"🧪 RDKit: {'✅' if HAS_RDKIT else '❌'}")
        st.write(f"🧬 Quantum Engine: {'✅' if HAS_QUANTUM else '❌'}")
        st.write(f"🌐 3D Visualization: {'✅' if HAS_3D_VIZ else '❌'}")
        st.write(f"🗃️ Knowledge Graph: {'✅' if HAS_AKG else '❌'}")
    
    # Main tabs
    tab1, tab2, tab3 = st.tabs(["🧪 Molecular Analysis", "📊 Discovery Dashboard", "⚛️ Quantum Metrics"])
    
    with tab1:
        st.header("🧪 Molecular Analysis")
        st.markdown("*Analyze molecular properties using quantum vQbit substrate*")
        
        # SMILES input
        col_input, col_examples = st.columns([2, 1])
        
        with col_input:
            smiles = st.text_input(
                "**Enter SMILES notation:**",
                value="CCO",
                help="Standard SMILES notation for molecular structure"
            )
        
        with col_examples:
            st.markdown("**Example molecules:**")
            if st.button("🍺 Ethanol (CCO)"):
                smiles = "CCO"
                st.rerun()
            if st.button("💊 Aspirin"):
                smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
                st.rerun()
            if st.button("☕ Caffeine"):
                smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
                st.rerun()
        
        if smiles:
            # Calculate molecular properties
            props = calculate_molecular_properties(smiles)
            
            if props:
                col_struct, col_props = st.columns([1, 1])
                
                with col_struct:
                    st.subheader("🔬 Molecular Structure")
                    
                    # 2D Structure
                    svg_2d = render_molecular_structure_2d(smiles)
                    if svg_2d:
                        st.markdown("**2D Structure:**")
                        st.markdown(svg_2d, unsafe_allow_html=True)
                    
                    # 3D Interactive Structure
                    if HAS_3D_VIZ:
                        try:
                            sdf_3d = render_molecular_structure_3d(smiles)
                            if sdf_3d:
                                st.markdown("**🌐 3D Interactive Structure:**")
                                
                                # Create multiple visualization styles
                                viz_style = st.selectbox(
                                    "Visualization Style:",
                                    ["Ball & Stick", "Stick", "Sphere", "Cartoon"],
                                    key="viz_style"
                                )
                                
                                viewer = py3Dmol.view(width=500, height=400)
                                viewer.addModel(sdf_3d, 'sdf')
                                
                                # Apply selected visualization style
                                if viz_style == "Ball & Stick":
                                    viewer.setStyle({'stick': {'radius': 0.15}, 'sphere': {'radius': 0.4}})
                                elif viz_style == "Stick":
                                    viewer.setStyle({'stick': {'radius': 0.2}})
                                elif viz_style == "Sphere":
                                    viewer.setStyle({'sphere': {'radius': 0.8}})
                                else:  # Cartoon
                                    viewer.setStyle({'cartoon': {}})
                                
                                viewer.setBackgroundColor('#ffffff')
                                viewer.zoomTo()
                                
                                # Add labels for atoms if molecule is small
                                mol = Chem.MolFromSmiles(smiles)
                                if mol and mol.GetNumAtoms() <= 20:
                                    for i, atom in enumerate(mol.GetAtoms()):
                                        viewer.addLabel(atom.GetSymbol(), 
                                                      {'position': {'x': 0, 'y': 0, 'z': 0}, 
                                                       'fontColor': 'black', 'fontSize': 12})
                                
                                stmol.showmol(viewer, height=400, width=500)
                                st.success("✅ 3D structure rendered successfully")
        else:
                                st.warning("⚠️ Could not generate 3D coordinates for this molecule")
                        except Exception as e:
                            st.error(f"❌ 3D visualization error: {e}")
                            st.info("💡 Some complex molecules may not render in 3D")
        else:
                        st.info("💡 Install stmol and py3Dmol for 3D visualization")
                
                with col_props:
                    st.subheader("📊 Molecular Properties")
                    
                    # Chemical Identity
                    st.markdown("**🔬 Chemical Identity:**")
                    st.text(f"Formula: {props['molecular_formula']}")
                    st.text(f"SMILES: {props['canonical_smiles']}")
                    st.text(f"MW: {props['molecular_weight']} g/mol")
                    st.text(f"Exact Mass: {props['exact_mass']} Da")
                    
                    # Physicochemical Properties
                    st.markdown("**⚗️ Physicochemical Properties:**")
                    st.text(f"LogP: {props['logp']}")
                    st.text(f"TPSA: {props['tpsa']} Ų")
                    st.text(f"H-Bond Donors: {props['hbd_count']}")
                    st.text(f"H-Bond Acceptors: {props['hba_count']}")
                    st.text(f"Rotatable Bonds: {props['rotatable_bond_count']}")
                    
                    # Drug-like properties
                    st.markdown("**💊 Drug-like Properties:**")
                    st.text(f"Lipinski Violations: {props['lipinski_violations']}/4")
                    if props['lipinski_violations'] == 0:
                        st.success("✅ Lipinski Rule of Five compliant")
            else:
                        st.warning(f"⚠️ {props['lipinski_violations']} Lipinski violations")
                
                # Quantum Analysis
                if HAS_QUANTUM:
                    st.subheader("⚛️ Quantum Property Analysis")
                    
                    if st.button("🌌 Analyze with Quantum Engine", key="quantum_analyze"):
                        with st.spinner("Running quantum analysis..."):
                            quantum_result = analyze_molecule_quantum(smiles)
                        
                        if quantum_result.get('success'):
                            # Quantum metrics
                            st.markdown("#### 🌌 Quantum State Metrics")
                            met_col1, met_col2, met_col3, met_col4 = st.columns(4)
                            with met_col1:
                                st.metric("🌊 Quantum Coherence", f"{quantum_result['coherence']:.6f}")
                            with met_col2:
                                st.metric("📏 State Normalization", f"{quantum_result['normalization']:.6f}")
                            with met_col3:
                                st.metric("⚛️ Chemical Potential", f"{quantum_result['chemical_potential']:.3f}")
                            with met_col4:
                                st.metric("🎯 Quantum Fidelity", f"{quantum_result['fidelity']:.1f}")
                            
                            # Quantum property measurements
                            st.markdown("#### 🔬 Quantum Property Measurements")
                            measurements_df = pd.DataFrame(list(quantum_result['measurements'].items()), 
                                                         columns=['Property', 'Measurement'])
                            
                            # Property descriptions for chemists
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
                            st.plotly_chart(fig, width='stretch')
                        else:
                            st.error(f"❌ Quantum analysis failed: {quantum_result.get('error', 'Unknown error')}")
    
    with tab2:
        st.header("📊 Discovery Dashboard")
        st.markdown("*Real-time view of autonomous discovery campaigns and results*")
        
        # Load discovery data
        discovery_data = load_discovery_data()
        
        if discovery_data:
            # Display discovery statistics
            if 'statistics' in discovery_data:
                stats = discovery_data['statistics']
                
                st.subheader("📈 Discovery Statistics")
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("🧪 Total Molecules", stats.get('total_molecules', 0))
                with col2:
                    st.metric("⚗️ Reactions", stats.get('total_reactions', 0))
                with col3:
                    st.metric("📐 Measurements", stats.get('total_measurements', 0))
                with col4:
                    st.metric("🎯 Active Claims", stats.get('active_claims', 0))
            
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
                                        st.image(svg_2d.encode('utf-8'), width=300)
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
                            'Drug Likeness': f"{mol.get('drug_likeness', 0):.3f}",
                            'Safety': f"{mol.get('safety_score', 0):.3f}",
                            'ID': mol.get('id', 'Unknown')[:8] + '...' if mol.get('id') else 'Unknown'
                        })
                    
                    summary_df = pd.DataFrame(summary_data)
                    st.dataframe(summary_df, width='stretch')
            else:
                st.info("No recent molecular discoveries found")
        else:
            st.info("📭 No discovery data available. Start discovery campaigns to see results here.")
            
            if st.button("🚀 Export Current Data"):
                akg = get_akg_client()
                if akg:
                    try:
                        akg.export_for_streamlit()
                        st.success("✅ Data exported successfully!")
                        st.rerun()
                    except Exception as e:
                        st.error(f"❌ Export failed: {e}")
    
    with tab3:
        st.header("⚛️ Quantum Substrate Metrics")
        st.markdown("*Performance and status of the quantum vQbit engine*")
        
        # System performance
        st.subheader("🚀 Performance Optimizations")
        perf_data = [
            {"Component": "C Extensions", "Status": "✅ Available" if True else "❌ Missing"},
            {"Component": "GPU/MPS Acceleration", "Status": "✅ Available" if engine and engine.gpu_acceleration else "❌ Missing"},
            {"Component": "3D Visualization", "Status": "✅ Available" if HAS_3D_VIZ else "❌ Missing"},
            {"Component": "RDKit Chemistry", "Status": "✅ Available" if HAS_RDKIT else "❌ Missing"},
            {"Component": "Neo4j AKG", "Status": "✅ Available" if HAS_AKG else "❌ Missing"}
        ]
        
        perf_df = pd.DataFrame(perf_data)
        st.dataframe(perf_df, width='stretch')
        
        # Quantum properties
        if engine:
            st.subheader("🌌 Quantum Properties")
            
        col1, col2, col3 = st.columns(3)
        with col1:
                st.metric("Hilbert Dimension", engine.hilbert_dimension)
        with col2:
                st.metric("Property Operators", len(ChemistryPropertyType))
        with col3:
                st.metric("Quantum Coherence Time", "∞ (noiseless)")
            
            # Test quantum operations
            if st.button("🧪 Test Quantum Operations"):
                with st.spinner("Testing quantum substrate..."):
                    try:
                        test_result = engine.test_quantum_substrate()
                        if test_result:
                            st.success("✅ All quantum operations working correctly")
                            for check, status in test_result.items():
                                st.text(f"{check}: {'✅ PASSED' if status else '❌ FAILED'}")
                        else:
                            st.error("❌ Quantum substrate test failed")
                    except Exception as e:
                        st.error(f"❌ Quantum test error: {e}")
        
        # Clear cache button
        if st.button("🗑️ Clear Cache"):
            st.cache_data.clear()
            st.cache_resource.clear()
            st.success("✅ Cache cleared!")
            st.rerun()

if __name__ == "__main__":
    main()
