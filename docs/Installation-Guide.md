# Installation Guide

> **🎯 NEW**: Complete problem-solution ontology system with **12,144 validated instances**!

Complete setup instructions for FoTChemistry quantum-guided molecular discovery system.

## 🎯 Prerequisites

### System Requirements

- **OS**: macOS 10.15+, Linux Ubuntu 18.04+, Windows 10+ (WSL2 recommended)
- **Python**: 3.9 or higher
- **Memory**: 8GB RAM minimum, 16GB+ recommended
- **Storage**: 5GB free space
- **GPU**: Apple Silicon (MPS) or CUDA-compatible GPU (optional but recommended)

### Required Dependencies

#### Neo4j Database
FoTChemistry requires a running Neo4j instance for the Agentic Knowledge Graph:

```bash
# Option 1: Neo4j Desktop (Recommended)
# Download from: https://neo4j.com/download/

# Option 2: Docker
docker run -d \
  --name neo4j-fotchem \
  -p 7474:7474 -p 7687:7687 \
  -e NEO4J_AUTH=neo4j/fotquantum \
  neo4j:latest

# Option 3: Homebrew (macOS)
brew install neo4j
neo4j start
```

**Important**: Set the Neo4j password to `fotquantum` for compatibility.

#### RDKit Chemistry Library
```bash
# Option 1: Conda (Recommended)
conda install -c conda-forge rdkit

# Option 2: pip (may require compilation)
pip install rdkit
```

## 📦 Installation Methods

### Method 1: Clone from GitHub (Recommended)

```bash
# Clone the repository
git clone https://github.com/FortressAI/FoTChemistry.git
cd FoTChemistry

# Create virtual environment
python -m venv fotchem-env
source fotchem-env/bin/activate  # On Windows: fotchem-env\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Compile quantum operations (optional, for performance)
cd core
python setup_quantum_ops.py build_ext --inplace
cd ..
```

### Method 2: pip install (Future)

```bash
# Will be available when published to PyPI
pip install fot-chemistry
```

## 🔧 Configuration

### 1. Neo4j Configuration

Ensure Neo4j is running with the correct credentials:

```bash
# Test Neo4j connection
python -c "
from neo4j import GraphDatabase
driver = GraphDatabase.driver('bolt://localhost:7687', auth=('neo4j', 'fotquantum'))
with driver.session() as session:
    result = session.run('RETURN 1 as test')
    print('✅ Neo4j connection successful')
driver.close()
"
```

### 2. GPU Acceleration (Optional)

#### Apple Silicon (MPS)
```bash
# Verify MPS availability
python -c "
import torch
if torch.backends.mps.is_available():
    print('✅ MPS acceleration available')
else:
    print('❌ MPS not available')
"
```

#### NVIDIA CUDA
```bash
# Verify CUDA availability
python -c "
import torch
if torch.cuda.is_available():
    print(f'✅ CUDA available: {torch.cuda.get_device_name()}')
else:
    print('❌ CUDA not available')
"
```

### 3. Verify Installation

Run the comprehensive verification script:

```bash
python -c "
import sys
print(f'Python: {sys.version}')

# Test core imports
try:
    from core.chemistry_vqbit_engine import ChemistryVQbitEngine
    print('✅ Quantum engine available')
except ImportError as e:
    print(f'❌ Quantum engine import failed: {e}')

try:
    from rdkit import Chem
    print('✅ RDKit available')
except ImportError:
    print('❌ RDKit not available')

try:
    import stmol, py3Dmol
    print('✅ 3D visualization available')
except ImportError:
    print('❌ 3D visualization not available')

try:
    from akg.client import AKG
    akg = AKG()
    print('✅ AKG connection successful')
except Exception as e:
    print(f'❌ AKG connection failed: {e}')
"
```

## 🚀 First Run

### 1. Start the Discovery Engine

```bash
# Test molecular discovery
python continuous_chemistry_discovery.py --test-mode --batch-size 3
```

Expected output:
```
✅ Chemistry Discovery Engine Initialized
✅ Quantum engine active: 8096-dimensional Hilbert space
✅ Real molecular generator initialized with RDKit
✅ Connected to EXISTING Neo4j instance
✅ Discovery 1: c1ccc(-c2cncnc2)cc1 (score: 0.746)
✅ Test completed: 3 discoveries
```

### 2. Launch Streamlit Dashboard

```bash
# Start the web interface
streamlit run streamlit_app.py --server.port 8505
```

Navigate to `http://localhost:8505` to access the dashboard.

### 3. Test Molecular Analysis

1. Go to the "🧪 Molecular Analysis" tab
2. Enter a SMILES string (e.g., `CCO` for ethanol)
3. Click "🌌 Analyze with Quantum Engine"
4. Verify quantum properties are displayed

## 🐛 Troubleshooting

### Common Issues

#### "ModuleNotFoundError: No module named 'rdkit'"
```bash
# Solution: Install RDKit via conda
conda install -c conda-forge rdkit
```

#### "AttributeError: 'ChemistryVQbitEngine' object has no attribute..."
```bash
# Solution: Ensure you're using the latest version
git pull origin main
pip install -r requirements.txt
```

#### "Neo.ClientError.Security.AuthenticationRateLimit"
```bash
# Solution: Verify Neo4j credentials
# Default: username=neo4j, password=fotquantum
```

#### "RuntimeWarning: divide by zero encountered in matmul"
This is expected for the quantum unitarity checks with large matrices. It doesn't affect functionality.

#### 3D Visualization not working
```bash
# Install missing dependencies
pip install stmol py3Dmol ipython_genutils
```

### Platform-Specific Issues

#### macOS Apple Silicon
- Install homebrew: `/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"`
- Install python: `brew install python@3.9`
- Install RDKit: `conda install -c conda-forge rdkit`

#### Linux
- Install system dependencies: `sudo apt-get install python3-dev build-essential`
- Install RDKit: `conda install -c conda-forge rdkit`

#### Windows
- Use WSL2 for best compatibility
- Install Anaconda for Windows
- Use conda for all chemistry dependencies

## 🔄 Updates

### Updating FoTChemistry

```bash
# Pull latest changes
git pull origin main

# Update dependencies
pip install -r requirements.txt --upgrade

# Recompile C extensions if needed
cd core
python setup_quantum_ops.py build_ext --inplace
cd ..
```

### Version Verification

```bash
python -c "
import json
with open('package.json', 'r') as f:
    version = json.load(f)['version']
print(f'FoTChemistry version: {version}')
"
```

## 🎯 Next Steps

Once installation is complete:

1. **[Quick Start Tutorial](Quick-Start-Tutorial)** - Learn the basics
2. **[Streamlit Dashboard Guide](Streamlit-Dashboard-Guide)** - Explore the interface
3. **[Molecular Discovery Pipeline](Molecular-Discovery-Pipeline)** - Run discovery campaigns
4. **[API Reference](API-Reference)** - Integrate with your code

## 📞 Support

If you encounter issues not covered here:

1. Check the [Troubleshooting](Troubleshooting) page
2. Search [GitHub Issues](https://github.com/FortressAI/FoTChemistry/issues)
3. Create a new issue with your system details and error messages

---

**🎉 Ready to discover molecules with quantum mechanics!**
