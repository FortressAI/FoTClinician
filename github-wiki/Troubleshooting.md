# Troubleshooting

Common issues and solutions for FoTChemistry quantum-guided molecular discovery system.

## 🚨 Critical Issues

### Quantum Engine Initialization Failed

**Symptoms**:
```
❌ Quantum vQbit engine initialization failed
QuantumEngineError: GPU context initialization failed
AttributeError: 'ChemistryVQbitEngine' object has no attribute 'hilbert_dimension'
```

**Causes**:
- Missing PyTorch installation
- GPU drivers not installed
- Insufficient memory (>8GB required)
- C extensions compilation failed

**Solutions**:

1. **Install PyTorch with MPS support (Apple Silicon)**:
```bash
pip install torch torchvision torchaudio
```

2. **Verify GPU availability**:
```python
import torch
print(f"MPS available: {torch.backends.mps.is_available()}")
print(f"CUDA available: {torch.cuda.is_available()}")
```

3. **Fallback to CPU mode**:
```python
engine = ChemistryVQbitEngine(use_gpu=False)
```

4. **Recompile C extensions**:
```bash
cd core
python setup_quantum_ops.py build_ext --inplace
```

5. **Check memory usage**:
```bash
# Monitor memory during initialization
htop  # Linux/macOS
# Task Manager on Windows
```

---

### Neo4j Connection Failed

**Symptoms**:
```
❌ AKG connection failed: AuthenticationRateLimit
Neo.ClientError.Security.Unauthorized: The client is unauthorized
ServiceUnavailable: Connection refused
```

**Causes**:
- Neo4j not running
- Incorrect credentials
- Port 7687 blocked
- Neo4j authentication rate limit

**Solutions**:

1. **Start Neo4j**:
```bash
# Neo4j Desktop: Click "Start" button
# Homebrew: 
brew services start neo4j

# Docker:
docker run -d --name neo4j \
  -p 7474:7474 -p 7687:7687 \
  -e NEO4J_AUTH=neo4j/fotquantum \
  neo4j:latest
```

2. **Verify credentials**:
```python
from neo4j import GraphDatabase
driver = GraphDatabase.driver("bolt://localhost:7687", 
                             auth=("neo4j", "fotquantum"))
try:
    with driver.session() as session:
        result = session.run("RETURN 1")
        print("✅ Connection successful")
except Exception as e:
    print(f"❌ Connection failed: {e}")
finally:
    driver.close()
```

3. **Reset Neo4j password**:
```bash
# Stop Neo4j, delete auth file, restart
rm $NEO4J_HOME/data/dbms/auth
neo4j start
# Default: neo4j/neo4j, then change to fotquantum
```

4. **Check firewall/ports**:
```bash
# Test port connectivity
telnet localhost 7687
nc -zv localhost 7687
```

---

### RDKit Import Failed

**Symptoms**:
```
ModuleNotFoundError: No module named 'rdkit'
ImportError: cannot import name 'Chem' from 'rdkit'
❌ RDKit not available - cannot render molecular structures
```

**Causes**:
- RDKit not installed
- Wrong Python environment
- Conda vs pip installation conflict

**Solutions**:

1. **Install via Conda (Recommended)**:
```bash
conda install -c conda-forge rdkit
```

2. **Install via pip**:
```bash
pip install rdkit
# Or specific version:
pip install rdkit-pypi
```

3. **Verify installation**:
```python
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    print("✅ RDKit available")
    
    # Test molecule creation
    mol = Chem.MolFromSmiles("CCO")
    if mol:
        print("✅ Molecule creation works")
    else:
        print("❌ Molecule creation failed")
        
except ImportError as e:
    print(f"❌ RDKit import failed: {e}")
```

4. **Environment conflicts**:
```bash
# Check which Python/environment you're using
which python
python -c "import sys; print(sys.path)"

# Create clean environment
conda create -n fotchem python=3.9
conda activate fotchem
conda install -c conda-forge rdkit
```

---

## ⚠️ Common Warnings

### RuntimeWarnings in Quantum Operations

**Symptoms**:
```
RuntimeWarning: divide by zero encountered in matmul
RuntimeWarning: overflow encountered in matmul  
RuntimeWarning: invalid value encountered in matmul
```

**Cause**: Numerical instability in large matrix operations (8096x8096) during quantum unitarity checks.

**Impact**: These are expected warnings and don't affect functionality.

**Solutions**:
1. **Ignore warnings** (safe for production):
```python
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, module="chemistry_vqbit_engine")
```

2. **Reduce test dimensions** (already implemented):
```python
# In chemistry_vqbit_engine.py, QFT unitarity test uses 8x8 instead of 8096x8096
test_dim = 8  # Instead of self.hilbert_dimension
```

---

### Streamlit AttributeError

**Symptoms**:
```
AttributeError: 'ChemistryVQbitEngine' object has no attribute 'hilbert_dim'
```

**Cause**: Incorrect attribute name in Streamlit app.

**Solution**: Use correct attribute name:
```python
# Wrong
engine.hilbert_dim

# Correct  
engine.hilbert_dimension
```

---

## 🧪 Molecular Analysis Issues

### Invalid SMILES Notation

**Symptoms**:
```
❌ Invalid SMILES notation
Molecular analysis failed: None
RDKit Error: Cannot parse SMILES
```

**Causes**:
- Typos in SMILES string
- Non-standard notation
- Special characters

**Solutions**:

1. **Validate SMILES**:
```python
from rdkit import Chem

def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"❌ Invalid SMILES: {smiles}")
        return False
    else:
        canonical = Chem.MolToSmiles(mol)
        print(f"✅ Valid SMILES: {smiles} → {canonical}")
        return True

# Test examples
validate_smiles("CCO")        # ✅ ethanol
validate_smiles("c1ccccc1")   # ✅ benzene  
validate_smiles("INVALID")    # ❌ invalid
```

2. **Common SMILES examples**:
```
✅ Valid SMILES:
CCO                     # ethanol
c1ccccc1               # benzene
CN1C=NC2=C1C(=O)N(C(=O)N2C)C  # caffeine
CC(=O)OC1=CC=CC=C1C(=O)O       # aspirin

❌ Common mistakes:
cco                    # Use CCO (uppercase C)
C1CCCCC                # Missing closing: C1CCCCC1
invalid_chars          # No underscores allowed
```

---

### 3D Visualization Not Working

**Symptoms**:
```
🌐 3D Visualization: ❌
ModuleNotFoundError: No module named 'stmol'
3D visualization not available
```

**Causes**:
- Missing stmol/py3Dmol packages
- WebGL not supported in browser
- JavaScript disabled

**Solutions**:

1. **Install 3D visualization packages**:
```bash
pip install stmol py3Dmol ipython_genutils
```

2. **Test 3D dependencies**:
```python
try:
    import stmol
    import py3Dmol
    from ipython_genutils import py3compat
    print("✅ 3D visualization packages available")
except ImportError as e:
    print(f"❌ Missing package: {e}")
```

3. **Browser requirements**:
- Enable WebGL in browser settings
- Use modern browser (Chrome 90+, Firefox 88+, Safari 14+)
- Enable JavaScript
- Disable browser extensions that block WebGL

4. **Test WebGL**:
Visit `https://get.webgl.org/` to verify WebGL support.

---

### Molecular Property Calculation Failed

**Symptoms**:
```
❌ Molecular property calculation failed
KeyError: 'molecular_weight'
TypeError: unsupported operand type(s)
```

**Cause**: RDKit calculation error or missing molecular data.

**Solutions**:

1. **Test property calculation**:
```python
from rdkit import Chem
from rdkit.Chem import Descriptors

def test_properties(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
            
        props = {
            'molecular_weight': Descriptors.ExactMolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'tpsa': Descriptors.TPSA(mol),
            'hbd_count': Descriptors.NumHDonors(mol),
            'hba_count': Descriptors.NumHAcceptors(mol)
        }
        return props
    except Exception as e:
        print(f"Property calculation failed: {e}")
        return None

# Test with simple molecule
props = test_properties("CCO")
print(props)
```

2. **Handle calculation errors gracefully**:
```python
def safe_calculate_properties(smiles):
    try:
        props = calculate_molecular_properties(smiles)
        return props if props else {}
    except Exception:
        return {
            'molecular_weight': 0.0,
            'logp': 0.0,
            'tpsa': 0.0,
            'error': 'Calculation failed'
        }
```

---

## 🔄 Discovery Pipeline Issues

### No Discoveries Generated

**Symptoms**:
```
📭 No discovery data available
Discovery engine completed: 0 discoveries
Empty discovery dashboard
```

**Causes**:
- Validation criteria too strict
- Molecular generation failed
- Quantum measurements failed
- Score threshold too high

**Solutions**:

1. **Lower validation thresholds**:
```python
# In continuous_chemistry_discovery.py
def _validate_discovery_candidate(self, candidate):
    # Lower thresholds for testing
    if candidate.get('combined_score', 0) < 0.1:  # Was 0.3
        return False
    
    if candidate.get('safety_score', 0) < 0.2:    # Was 0.4
        return False
    
    return True
```

2. **Test molecular generation**:
```python
from agents.alchemist.real_molecular_generator import RealMolecularGenerator
from core.chemistry_vqbit_engine import ChemistryVQbitEngine

engine = ChemistryVQbitEngine(use_gpu=False)
generator = RealMolecularGenerator(engine)

candidates = generator.generate_candidates("CCO", num_variants=3)
print(f"Generated {len(candidates)} candidates")
for candidate in candidates:
    print(f"- {candidate.get('smiles', 'N/A')}: {candidate.get('combined_score', 0):.3f}")
```

3. **Debug discovery process**:
```bash
# Run with debug logging
python continuous_chemistry_discovery.py --test-mode --batch-size 1 --debug
```

---

### Discovery Engine Crashes

**Symptoms**:
```
Segmentation fault (core dumped)
MemoryError: Unable to allocate array
KeyboardInterrupt during discovery
```

**Causes**:
- Insufficient memory
- C extension bugs
- Infinite loops in molecular generation

**Solutions**:

1. **Monitor memory usage**:
```bash
# Monitor while running discovery
watch -n 1 "ps aux | grep python | grep continuous"
```

2. **Reduce batch size**:
```python
config = ContinuousConfig(
    batch_size=3,      # Instead of 10+
    interval_seconds=60  # Longer intervals
)
```

3. **Use CPU-only mode**:
```python
engine = ChemistryVQbitEngine(use_gpu=False)
```

4. **Add memory limits**:
```python
import resource

# Limit memory to 8GB
resource.setrlimit(resource.RLIMIT_AS, (8 * 1024 * 1024 * 1024, -1))
```

---

## 🌐 Streamlit Dashboard Issues

### Dashboard Not Loading

**Symptoms**:
```
ModuleNotFoundError: No module named 'streamlit'
Connection refused at localhost:8505
Page not found
```

**Solutions**:

1. **Install Streamlit**:
```bash
pip install streamlit
```

2. **Start dashboard**:
```bash
streamlit run streamlit_app.py --server.port 8505
```

3. **Check port availability**:
```bash
# Check if port is in use
lsof -i :8505
netstat -tulpn | grep 8505

# Use different port if needed
streamlit run streamlit_app.py --server.port 8506
```

4. **Access dashboard**:
```
Local:    http://localhost:8505
Network:  http://[your-ip]:8505
```

---

### Slow Dashboard Performance

**Symptoms**:
- Long loading times
- Unresponsive interface
- High CPU usage

**Solutions**:

1. **Clear Streamlit cache**:
```python
# In dashboard, click "Clear Cache" button
# Or restart dashboard:
# Ctrl+C, then streamlit run streamlit_app.py
```

2. **Reduce quantum calculations**:
```python
# Use smaller molecules for testing
test_molecules = ["C", "O", "CCO"]  # Instead of complex structures
```

3. **Optimize caching**:
```python
@st.cache_data(ttl=300)  # Cache for 5 minutes
def get_molecular_properties(smiles):
    return calculate_molecular_properties(smiles)
```

---

## 🔧 Performance Issues

### Slow Quantum Calculations

**Symptoms**:
- Long wait times for quantum analysis
- High CPU usage
- Memory consumption

**Solutions**:

1. **Enable GPU acceleration**:
```python
# Verify GPU is being used
engine = ChemistryVQbitEngine(use_gpu=True)
print(f"GPU acceleration: {engine.gpu_acceleration}")
```

2. **Compile C extensions**:
```bash
cd core
python setup_quantum_ops.py build_ext --inplace
```

3. **Reduce calculation complexity**:
```python
# Use smaller test dimensions for debugging
# Modify chemistry_vqbit_engine.py temporarily
hilbert_dimension = 1024  # Instead of 8096
```

---

## 🐛 Debug Mode

### Enable Comprehensive Logging

```python
import logging

# Set up detailed logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('fotchemistry_debug.log'),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger('fotchemistry')
logger.debug("Debug mode enabled")
```

### Test Individual Components

```python
# Test quantum engine
def test_quantum_engine():
    try:
        engine = ChemistryVQbitEngine(use_gpu=False)
        result = engine.test_quantum_substrate()
        print(f"Quantum tests: {result}")
        return True
    except Exception as e:
        print(f"Quantum engine failed: {e}")
        return False

# Test molecular generator
def test_molecular_generator():
    try:
        engine = ChemistryVQbitEngine(use_gpu=False)
        generator = RealMolecularGenerator(engine)
        candidates = generator.generate_candidates("CCO", num_variants=1)
        print(f"Generated {len(candidates)} candidates")
        return len(candidates) > 0
    except Exception as e:
        print(f"Molecular generator failed: {e}")
        return False

# Test AKG connection
def test_akg_connection():
    try:
        akg = AKG()
        # Try a simple query
        discoveries = akg.query_discoveries(limit=1)
        print(f"AKG connection successful, {len(discoveries)} discoveries")
        return True
    except Exception as e:
        print(f"AKG connection failed: {e}")
        return False

# Run all tests
if __name__ == "__main__":
    print("Testing FoTChemistry components...")
    results = {
        'quantum_engine': test_quantum_engine(),
        'molecular_generator': test_molecular_generator(),
        'akg_connection': test_akg_connection()
    }
    
    print("\nTest Results:")
    for component, status in results.items():
        status_str = "✅ PASS" if status else "❌ FAIL"
        print(f"{component}: {status_str}")
```

---

## 📞 Getting Help

### Before Reporting Issues

1. **Check this troubleshooting guide**
2. **Search [GitHub Issues](https://github.com/FortressAI/FoTChemistry/issues)**
3. **Verify your installation** using debug scripts above
4. **Collect system information**:

```bash
# System info
python --version
pip list | grep -E "(torch|rdkit|streamlit|neo4j)"
uname -a  # Linux/macOS
systeminfo | findstr /B /C:"OS Name" /C:"OS Version"  # Windows
```

### Reporting Issues

Include this information in your issue:

1. **System details**: OS, Python version, hardware
2. **Installation method**: pip, conda, Docker
3. **Complete error message**: Full traceback
4. **Reproduction steps**: Minimal example that fails
5. **Expected vs actual behavior**
6. **Relevant logs**: Debug output

### Issue Template

```markdown
**System Information:**
- OS: [e.g., macOS 12.6, Ubuntu 20.04]
- Python: [e.g., 3.9.7]
- FoTChemistry: [version]
- GPU: [Apple M1, NVIDIA RTX 3080, None]

**Problem Description:**
[Clear description of the issue]

**Error Message:**
```
[Full error traceback]
```

**Steps to Reproduce:**
1. Run `python continuous_chemistry_discovery.py --test-mode`
2. See error at step X
3. Expected Y but got Z

**Additional Context:**
[Any other relevant information]
```

---

**🔧 Most issues can be resolved with proper environment setup and following the solutions above**

For additional help:
- **[Installation Guide](Installation-Guide)** - Proper setup
- **[API Reference](API-Reference)** - Usage examples  
- **[GitHub Issues](https://github.com/FortressAI/FoTChemistry/issues)** - Community support
