# Streamlit Dashboard Guide

> **📊 MASSIVE DATASET**: Explore **6,443 real molecular discoveries** through interactive visualization!

Complete guide to using FoTChemistry's interactive web dashboard for molecular analysis, discovery monitoring, and quantum property exploration.

## 🌐 Dashboard Overview

The FoTChemistry Streamlit dashboard provides a comprehensive interface for:
- **🧪 Molecular Analysis**: Interactive quantum property analysis with 3D visualization
- **📊 Discovery Dashboard**: Real-time monitoring of autonomous molecular discoveries  
- **⚛️ Quantum Metrics**: Performance monitoring and system status

**Access**: `http://localhost:8505` (after running `streamlit run streamlit_app.py --server.port 8505`)

## 🧪 Molecular Analysis Tab

### Input Methods

#### 1. Manual SMILES Entry
```
Enter any valid SMILES notation:
- CCO (ethanol)
- c1ccccc1 (benzene)  
- CN1C=NC2=C1C(=O)N(C(=O)N2C)C (caffeine)
```

#### 2. Quick Examples
Click any example button for instant analysis:
- **🍺 Ethanol (CCO)**: Simple alcohol
- **💊 Aspirin**: Classic pharmaceutical
- **☕ Caffeine**: Stimulant alkaloid

### Molecular Structure Visualization

#### 2D Structure Display
- **Automatic Generation**: Creates clean 2D chemical diagrams
- **SVG Format**: Scalable vector graphics for crisp display
- **Standard Layout**: Following chemical drawing conventions

#### 3D Interactive Visualization

**Available Styles**:
- **Ball & Stick**: Shows atomic spheres connected by bonds (best for molecular geometry)
- **Stick**: Clean bond-only representation (good for complex molecules)
- **Sphere**: Space-filling atomic representation (shows molecular size)
- **Cartoon**: Simplified schematic view

**Interactive Controls**:
- **Rotate**: Click and drag to rotate molecule
- **Zoom**: Mouse wheel to zoom in/out
- **Labels**: Automatic atom labeling for molecules ≤20 atoms

**Technical Details**:
- **Engine**: py3Dmol + stmol integration
- **3D Coordinates**: RDKit UFF optimization
- **File Format**: SDF molecular data
- **Rendering**: WebGL hardware acceleration

### Molecular Properties

#### Chemical Identity
```
Formula: C8H9NO2 (aspirin example)
SMILES: CC(=O)OC1=CC=CC=C1C(=O)O
MW: 180.16 g/mol
Exact Mass: 180.0422 Da
```

#### Physicochemical Properties
```
LogP: 1.25 (lipophilicity)
TPSA: 63.6 Ų (polar surface area)
H-Bond Donors: 1
H-Bond Acceptors: 4
Rotatable Bonds: 3
```

#### Drug-like Properties
```
Lipinski Violations: 0/4
✅ Lipinski Rule of Five compliant
- MW ≤ 500 Da ✓
- LogP ≤ 5 ✓  
- H-Donors ≤ 5 ✓
- H-Acceptors ≤ 10 ✓
```

### Quantum Property Analysis

#### Quantum State Creation
Click **"🌌 Analyze with Quantum Engine"** to:
1. Create 8096-dimensional quantum vQbit state
2. Apply chemistry-specific property operators
3. Measure quantum expectation values
4. Calculate coherence and fidelity metrics

#### Quantum Metrics Display

**State Metrics**:
- **🌊 Quantum Coherence**: ~0.000199 (L1 coherence measure)
- **📏 State Normalization**: 1.000000 (quantum state normalization)
- **⚛️ Chemical Potential**: Variable (molecular energy landscape)
- **🎯 Quantum Fidelity**: 0.2-1.0 (state quality measure)

**Property Measurements**:
- **Bioactivity**: 0.0-1.0 (biological activity potential)
- **Sustainability**: 0.0-1.0 (environmental impact score)
- **Reproducibility**: 0.0-1.0 (experimental reliability)
- **Efficiency**: 0.0-1.0 (synthetic efficiency score)

#### Interactive Visualization
- **Bar Chart**: Property measurements with color coding
- **Hover Details**: Property descriptions and interpretations
- **Plasma Colormap**: Visual gradient from low to high values
- **Real-time Updates**: Immediate results from quantum engine

## 📊 Discovery Dashboard Tab

### Discovery Statistics

**Metrics Overview**:
```
🧪 Total Molecules: 5 (recent discoveries)
⚗️ Reactions: 15 (transformations attempted)  
📐 Measurements: 20 (quantum measurements)
🎯 Active Claims: 3 (pending validations)
```

### Recent Molecular Discoveries

#### Discovery Cards
Each discovery shows:
- **Molecular Structure**: 2D chemical diagram
- **SMILES Notation**: Chemical identifier
- **Combined Score**: Overall molecular quality (0.0-1.0)
- **Drug Likeness**: Pharmaceutical potential (0.0-1.0)
- **Safety Score**: Predicted safety profile (0.0-1.0)  
- **Quantum Coherence**: Quantum state stability
- **Discovery ID**: Unique identifier for tracking

#### Example Discovery Entry
```
🧬 Discovery 1: c1ccc(-c2cncnc2)cc1 (Score: 0.746)

📋 Discovery Details:
SMILES: c1ccc(-c2cncnc2)cc1
Combined Score: 0.746
Drug Likeness: 0.820
Safety Score: 0.750
Quantum Coherence: 0.000199
Discovery ID: chem_disc...

📊 Molecular Properties:
Formula: C11H8N2
MW: 168.2 g/mol
LogP: 2.85
```

### Discovery Summary Table

For campaigns with >5 discoveries:
```
SMILES                    | Score | Drug Likeness | Safety | ID
c1ccc(-c2cncnc2)cc1      | 0.746 | 0.820        | 0.750  | chem_disc...
CC(CN)NCC(O)c1ccc(O)...  | 0.724 | 0.780        | 0.860  | chem_disc...
c1ccc(C2CCCC2)cc1        | 0.756 | 0.710        | 0.920  | chem_disc...
```

### Data Export & Refresh

**Export Current Data**:
- Click **"🚀 Export Current Data"** to update the discovery file
- Exports Neo4j discoveries to `results/chemistry_discoveries.json`
- Enables cloud deployment data synchronization

**Automatic Refresh**:
- Dashboard automatically loads latest discovery data
- Local deployment: Real-time Neo4j queries
- Cloud deployment: Periodic JSON file updates

## ⚛️ Quantum Metrics Tab

### Performance Optimization Status

```
Component                | Status
C Extensions            | ✅ Available
GPU/MPS Acceleration    | ✅ Available
3D Visualization        | ✅ Available
RDKit Chemistry         | ✅ Available
Neo4j AKG              | ✅ Available
```

### Quantum Properties

**Core Metrics**:
- **Hilbert Dimension**: 8096 (quantum state space size)
- **Property Operators**: 4 (Bioactivity, Sustainability, Reproducibility, Efficiency)
- **Quantum Coherence Time**: ∞ (noiseless substrate)

### System Testing

**Quantum Operations Test**:
Click **"🧪 Test Quantum Operations"** to verify:
- Hilbert dimension integrity (8096 confirmed)
- Coherence time validation (infinite confirmed)
- Quantum fidelity verification (perfect confirmed)
- Decoherence rate check (zero confirmed)
- Identity operator validation (unitary confirmed)

**Expected Results**:
```
✅ All quantum operations working correctly
Hilbert dimension: ✅ PASSED
Coherence time: ✅ PASSED
Quantum fidelity: ✅ PASSED
Decoherence rate: ✅ PASSED
Identity operator: ✅ PASSED
```

### Cache Management

**Clear Cache**: 
- **Data Cache**: Clears molecular property calculations
- **Resource Cache**: Reloads quantum engine and AKG connections
- **Use Case**: After system updates or configuration changes

## 🔧 Sidebar System Status

### Quantum vQbit Engine
```
✅ Quantum Substrate Active
🌌 Hilbert Dimension: 8096
🚀 GPU Acceleration: Yes
```

### Agentic Knowledge Graph
```
✅ AKG Connected (Neo4j running)
❌ AKG Offline (connection failed)
```

### Feature Availability
```
🧪 RDKit: ✅
🧬 Quantum Engine: ✅
🌐 3D Visualization: ✅
🗃️ Knowledge Graph: ✅
```

## 🎨 UI/UX Features

### Responsive Design
- **Wide Layout**: Optimized for desktop/laptop screens
- **Column Layout**: Efficient space utilization
- **Mobile Friendly**: Responsive components for tablets

### Interactive Elements
- **Expandable Sections**: Click to expand discovery details
- **Hover Tooltips**: Property descriptions and help text
- **Progress Spinners**: Visual feedback during quantum calculations
- **Color Coding**: Status indicators and metric visualization

### Performance Features
- **Caching**: `@st.cache_data` and `@st.cache_resource` decorators
- **Lazy Loading**: Components load on-demand
- **Error Handling**: Graceful degradation for missing features
- **Background Processing**: Non-blocking quantum calculations

## 🚀 Advanced Usage

### Custom SMILES Collections

Create your own molecular libraries:
```python
pharmaceutical_molecules = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen  
    "CC(=O)NC1=CC=C(C=C1)O",     # Paracetamol
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
]

for smiles in pharmaceutical_molecules:
    # Analyze each molecule systematically
    analyze_molecule(smiles)
```

### Batch Analysis Scripts

```python
# Analyze multiple molecules programmatically
import streamlit as st
from streamlit_app import analyze_molecule_quantum

molecules = ["CCO", "c1ccccc1", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
results = []

for smiles in molecules:
    result = analyze_molecule_quantum(smiles)
    results.append({
        'smiles': smiles,
        'quantum_result': result
    })

# Display results in custom format
for result in results:
    st.write(f"Molecule: {result['smiles']}")
    st.write(f"Quantum Properties: {result['quantum_result']}")
```

### Integration with Discovery Engine

```python
# Monitor discoveries in real-time
from continuous_chemistry_discovery import ContinuousChemistryDiscoveryEngine

# Run discovery in background
engine = ContinuousChemistryDiscoveryEngine()
discoveries = engine.run_discovery_batch()

# Display in dashboard
for discovery in discoveries:
    st.write(f"New Discovery: {discovery.smiles}")
    # Automatic addition to dashboard data
```

## 🐛 Troubleshooting

### Common Issues

#### "Quantum analysis failed"
**Symptoms**: Error message during quantum analysis
**Causes**: 
- Invalid SMILES notation
- Quantum engine initialization failure
- Memory limitations

**Solutions**:
```python
# Verify SMILES validity
from rdkit import Chem
mol = Chem.MolFromSmiles("YOUR_SMILES_HERE")
if mol is None:
    print("Invalid SMILES")

# Check quantum engine
engine = get_quantum_engine()
if engine is None:
    print("Quantum engine failed to initialize")
```

#### "3D visualization not available"
**Symptoms**: Missing 3D molecular viewer
**Causes**: Missing stmol/py3Dmol packages
**Solutions**:
```bash
pip install stmol py3Dmol ipython_genutils
```

#### "No discoveries found"  
**Symptoms**: Empty discovery dashboard
**Causes**: Discovery engine not run or data not exported
**Solutions**:
```bash
# Run discovery engine first
python continuous_chemistry_discovery.py --test-mode --batch-size 3

# Then refresh dashboard
```

#### Slow Performance
**Symptoms**: Long loading times for quantum analysis
**Causes**: No GPU acceleration or large molecules
**Solutions**:
- Enable MPS/CUDA acceleration
- Try smaller molecules first
- Check system resources

### Debug Mode

Enable detailed logging:
```python
import logging
logging.basicConfig(level=logging.DEBUG)

# Run dashboard with verbose output
streamlit run streamlit_app.py --server.port 8505 --logger.level debug
```

### Browser Compatibility

**Recommended Browsers**:
- Chrome 90+ ✅
- Firefox 88+ ✅  
- Safari 14+ ✅
- Edge 90+ ✅

**WebGL Requirements**:
- Hardware acceleration enabled
- WebGL 2.0 support for 3D visualization

## 📱 Cloud Deployment

### Streamlit Cloud Setup

1. **Repository**: Push to GitHub with all dependencies
2. **Requirements**: Ensure `streamlit_requirements.txt` includes all packages
3. **Data**: Export discoveries to JSON for cloud access
4. **Configuration**: Set environment variables if needed

### Local vs Cloud Differences

| Feature | Local | Cloud |
|---------|-------|-------|
| **Neo4j Access** | ✅ Direct connection | ❌ Not available |
| **Discovery Data** | Real-time from Neo4j | JSON file export |
| **3D Visualization** | ✅ Full WebGL | ✅ Full WebGL |
| **Quantum Engine** | ✅ Full functionality | ⚠️ Limited (no C extensions) |
| **Performance** | Optimal | Good |

## 📚 Related Documentation

- **[Quick Start Tutorial](Quick-Start-Tutorial)** - Get started in 10 minutes
- **[Molecular Discovery Pipeline](Molecular-Discovery-Pipeline)** - Understanding discoveries
- **[API Reference](API-Reference)** - Integration and automation
- **[Troubleshooting](Troubleshooting)** - Detailed problem resolution

---

**🌐 Interactive exploration of quantum-guided molecular discovery**
