# API Reference

Complete API documentation for FoTChemistry's quantum-guided molecular discovery system.

## 🏗️ Core Architecture

### Main Components

```python
from core.chemistry_vqbit_engine import ChemistryVQbitEngine, ChemistryPropertyType
from agents.alchemist.real_molecular_generator import RealMolecularGenerator
from continuous_chemistry_discovery import ContinuousChemistryDiscoveryEngine
from akg.client import AKG
```

## ⚛️ ChemistryVQbitEngine

### Class Definition

```python
class ChemistryVQbitEngine:
    """8096-dimensional quantum vQbit engine for molecular property analysis"""
    
    def __init__(self, use_gpu: bool = True):
        """Initialize quantum engine with optional GPU acceleration"""
        
    @property
    def hilbert_dimension(self) -> int:
        """Get Hilbert space dimension (8096)"""
        
    @property 
    def gpu_acceleration(self) -> bool:
        """Check if GPU acceleration is active"""
```

### Core Methods

#### create_molecular_vqbit()
```python
def create_molecular_vqbit(self, property_scores: Dict[ChemistryPropertyType, float]) -> MolecularVQbitState:
    """
    Create quantum vQbit state from molecular property scores
    
    Args:
        property_scores: Dict mapping property types to scores [0.0, 1.0]
        
    Returns:
        MolecularVQbitState: Quantum state with amplitudes and metadata
        
    Example:
        scores = {
            ChemistryPropertyType.BIOACTIVITY: 0.8,
            ChemistryPropertyType.SUSTAINABILITY: 0.6,
            ChemistryPropertyType.REPRODUCIBILITY: 0.9,
            ChemistryPropertyType.EFFICIENCY: 0.7
        }
        vqbit_state = engine.create_molecular_vqbit(scores)
    """
```

#### measure_property()
```python
def measure_property(self, vqbit_state: MolecularVQbitState, 
                    property_type: ChemistryPropertyType) -> float:
    """
    Measure quantum property expectation value
    
    Args:
        vqbit_state: Quantum state to measure
        property_type: Which property to measure
        
    Returns:
        float: Measured expectation value [0.0, 1.0]
        
    Example:
        bioactivity = engine.measure_property(vqbit_state, ChemistryPropertyType.BIOACTIVITY)
    """
```

#### test_quantum_substrate()
```python
def test_quantum_substrate(self) -> Dict[str, bool]:
    """
    Verify quantum substrate integrity
    
    Returns:
        Dict[str, bool]: Test results for each quantum property
        
    Example:
        results = engine.test_quantum_substrate()
        # {'Hilbert dimension': True, 'Coherence time': True, ...}
    """
```

### Data Structures

#### ChemistryPropertyType
```python
class ChemistryPropertyType(Enum):
    """Quantum property types for molecular analysis"""
    BIOACTIVITY = "bioactivity"          # Biological activity potential
    SUSTAINABILITY = "sustainability"    # Environmental impact
    REPRODUCIBILITY = "reproducibility"  # Experimental reliability  
    EFFICIENCY = "efficiency"            # Synthetic efficiency
```

#### MolecularVQbitState
```python
@dataclass
class MolecularVQbitState:
    """Quantum state representation for molecules"""
    amplitudes: np.ndarray              # Complex amplitudes (8096-dim)
    property_scores: Dict               # Initial property scores
    hilbert_dimension: int = 8096       # State space dimension
    coherence: Optional[float] = None   # L1 coherence measure
    fidelity: Optional[float] = None    # State fidelity
```

## 🧬 RealMolecularGenerator

### Class Definition

```python
class RealMolecularGenerator:
    """Real molecular generation using RDKit and quantum guidance"""
    
    def __init__(self, quantum_engine: ChemistryVQbitEngine, 
                 config: Optional[Dict[str, Any]] = None):
        """Initialize with quantum engine and optional configuration"""
```

### Core Methods

#### generate_candidates()
```python
def generate_candidates(self, seed_smiles: str, num_variants: int = 5, 
                       objective: str = "drug_discovery") -> List[Dict[str, Any]]:
    """
    Generate molecular candidates using quantum guidance
    
    Args:
        seed_smiles: Starting molecular structure (SMILES)
        num_variants: Number of variants to generate
        objective: Discovery objective ('drug_discovery', 'green_chemistry')
        
    Returns:
        List[Dict]: Generated molecular candidates with scores
        
    Example:
        generator = RealMolecularGenerator(quantum_engine)
        candidates = generator.generate_candidates("CCO", num_variants=5)
        
        # Each candidate contains:
        # {
        #     'smiles': 'generated_smiles',
        #     'combined_score': 0.75,
        #     'drug_likeness': 0.82,
        #     'safety_score': 0.78,
        #     'quantum_coherence': 0.000199,
        #     'modifications': [...],
        #     'is_valid': True
        # }
    """
```

### Generation Strategies

#### Structural Modifications
```python
def _generate_structural_modifications(self, mol: Chem.Mol, count: int) -> List[Dict]:
    """Apply chemical transformations (functional group changes, etc.)"""
```

#### Bioisosteric Replacements  
```python
def _generate_bioisosteric_replacements(self, mol: Chem.Mol, count: int) -> List[Dict]:
    """Replace functional groups with bioisosteres"""
```

#### Fragment Additions
```python
def _generate_fragment_additions(self, mol: Chem.Mol, count: int) -> List[Dict]:
    """Add molecular fragments to create new structures"""
```

## 🔄 ContinuousChemistryDiscoveryEngine

### Class Definition

```python
class ContinuousChemistryDiscoveryEngine:
    """Autonomous molecular discovery engine"""
    
    def __init__(self, config: ContinuousConfig):
        """Initialize with discovery configuration"""
```

### Configuration

#### ContinuousConfig
```python
@dataclass
class ContinuousConfig:
    """Configuration for continuous discovery"""
    batch_size: int = 10                    # Molecules per batch
    interval_seconds: int = 30              # Time between batches
    objectives: List[str] = ["drug_discovery"]  # Discovery objectives
    min_score_threshold: float = 0.3        # Minimum quality threshold
    max_memory_gb: float = 16.0             # Memory limit
    max_cpu_percent: float = 85.0           # CPU usage limit
    output_directory: str = "results"       # Output location
    seed_molecules: List[str] = field(      # Starting molecules
        default_factory=lambda: [
            "CCO", "c1ccccc1", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "CC(C)NCC(C1=CC(=C(C=C1)O)CO)O", "CC(=O)OC1=CC=CC=C1C(=O)O"
        ]
    )
```

### Core Methods

#### run_discovery_batch()
```python
def run_discovery_batch(self, batch_id: Optional[str] = None) -> List[ChemicalDiscovery]:
    """
    Run single discovery batch
    
    Args:
        batch_id: Optional batch identifier
        
    Returns:
        List[ChemicalDiscovery]: Discovered molecules
        
    Example:
        engine = ContinuousChemistryDiscoveryEngine(config)
        discoveries = engine.run_discovery_batch("batch_001")
    """
```

#### run_continuous_discovery()
```python
def run_continuous_discovery(self) -> None:
    """
    Run continuous discovery loop until stopped
    
    Example:
        engine = ContinuousChemistryDiscoveryEngine(config)
        engine.run_continuous_discovery()  # Runs indefinitely
    """
```

### Data Structures

#### ChemicalDiscovery
```python
@dataclass  
class ChemicalDiscovery:
    """Represents a discovered molecule"""
    id: str                          # Unique identifier
    smiles: str                      # SMILES notation
    combined_score: float            # Overall quality [0.0, 1.0]
    drug_likeness: float            # Drug-like properties [0.0, 1.0]
    safety_score: float             # Safety assessment [0.0, 1.0]
    quantum_coherence: float        # Quantum coherence measure
    molecular_weight: float         # Molecular weight (Da)
    logp: float                     # Lipophilicity
    tpsa: float                     # Topological polar surface area
    lipinski_violations: int        # Rule of Five violations
    timestamp: datetime             # Discovery time
    batch_id: str                   # Batch identifier
    objective: str                  # Discovery objective
    generation_strategy: str        # How it was generated
    validation_passed: bool         # Validation status
```

## 🗃️ AKG (Agentic Knowledge Graph)

### Class Definition

```python
class AKG:
    """Neo4j-based Agentic Knowledge Graph client"""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize with Neo4j connection configuration"""
```

### Core Methods

#### store_discovery()
```python
def store_discovery(self, discovery: ChemicalDiscovery) -> bool:
    """
    Store discovery in Neo4j knowledge graph
    
    Args:
        discovery: ChemicalDiscovery object to store
        
    Returns:
        bool: Success status
        
    Example:
        akg = AKG()
        success = akg.store_discovery(discovery)
    """
```

#### query_discoveries()
```python
def query_discoveries(self, limit: int = 50, 
                     objective: Optional[str] = None) -> List[Dict]:
    """
    Query recent discoveries from knowledge graph
    
    Args:
        limit: Maximum number of discoveries to return
        objective: Filter by discovery objective
        
    Returns:
        List[Dict]: Discovery data
        
    Example:
        akg = AKG()
        recent = akg.query_discoveries(limit=10, objective="drug_discovery")
    """
```

#### export_for_streamlit()
```python  
def export_for_streamlit(self, output_file: str = "results/chemistry_discoveries.json") -> None:
    """
    Export discoveries to JSON for Streamlit dashboard
    
    Args:
        output_file: Path to output JSON file
        
    Example:
        akg = AKG()
        akg.export_for_streamlit("results/latest_discoveries.json")
    """
```

### Neo4j Schema

#### Node Types
```cypher
// Discovery nodes
(:FoTChem_Discovery {
    id: STRING,
    smiles: STRING,
    combined_score: FLOAT,
    drug_likeness: FLOAT,
    safety_score: FLOAT,
    quantum_coherence: FLOAT,
    timestamp: DATETIME,
    batch_id: STRING,
    objective: STRING
})

// Molecule nodes  
(:FoTChem_Molecule {
    id: STRING,
    smiles: STRING,
    molecular_formula: STRING,
    molecular_weight: FLOAT,
    inchi_key: STRING
})

// Campaign nodes
(:FoTChem_Campaign {
    id: STRING,
    name: STRING,
    objectives: LIST<STRING>,
    start_time: DATETIME,
    status: STRING
})
```

#### Relationships
```cypher
// Discovery relationships
(:FoTChem_Discovery)-[:DISCOVERED_MOLECULE]->(:FoTChem_Molecule)
(:FoTChem_Discovery)-[:PART_OF_CAMPAIGN]->(:FoTChem_Campaign)
(:FoTChem_Discovery)-[:GENERATED_FROM]->(:FoTChem_Molecule)  // seed molecule
```

## 🧪 Molecular Analysis Functions

### RDKit Integration

#### calculate_molecular_properties()
```python
def calculate_molecular_properties(smiles: str) -> Dict[str, Any]:
    """
    Calculate comprehensive molecular properties
    
    Args:
        smiles: SMILES notation string
        
    Returns:
        Dict with molecular properties
        
    Example:
        props = calculate_molecular_properties("CCO")
        # {
        #     'molecular_formula': 'C2H6O',
        #     'molecular_weight': 46.07,
        #     'logp': -0.31,
        #     'tpsa': 20.23,
        #     'hbd_count': 1,
        #     'hba_count': 1,
        #     'lipinski_violations': 0,
        #     'canonical_smiles': 'CCO',
        #     'inchi_key': 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
        # }
    """
```

#### render_molecular_structure_2d()
```python
def render_molecular_structure_2d(smiles: str) -> Optional[str]:
    """
    Generate 2D molecular structure as SVG
    
    Args:
        smiles: SMILES notation string
        
    Returns:
        Optional[str]: SVG string or None if failed
        
    Example:
        svg = render_molecular_structure_2d("CCO")
        # Returns SVG string for 2D molecule display
    """
```

#### render_molecular_structure_3d()
```python
def render_molecular_structure_3d(smiles: str) -> Optional[str]:
    """
    Generate 3D molecular coordinates as SDF
    
    Args:
        smiles: SMILES notation string
        
    Returns:
        Optional[str]: SDF string or None if failed
        
    Example:
        sdf = render_molecular_structure_3d("CCO")
        # Returns SDF string for 3D visualization
    """
```

## 🚀 Command Line Interface

### continuous_chemistry_discovery.py

```bash
# Basic usage
python continuous_chemistry_discovery.py [OPTIONS]

# Options:
--test-mode              # Run single batch for testing
--batch-size INT         # Molecules per batch (default: 10)
--interval INT           # Seconds between batches (default: 30)
--objectives LIST        # Discovery objectives (default: ["drug_discovery"])
--min-score FLOAT        # Minimum quality threshold (default: 0.3)
--continuous            # Run continuously until stopped
--output-dir PATH       # Output directory (default: "results")

# Examples:
python continuous_chemistry_discovery.py --test-mode --batch-size 5
python continuous_chemistry_discovery.py --continuous --interval 60 --objectives drug_discovery green_chemistry
```

### streamlit_app.py

```bash
# Launch dashboard
streamlit run streamlit_app.py [OPTIONS]

# Options:
--server.port INT       # Port number (default: 8501)
--server.address STR    # Server address (default: localhost)

# Example:
streamlit run streamlit_app.py --server.port 8505
```

## 🔧 Configuration Files

### Discovery Configuration (YAML)

```yaml
# discovery_config.yaml
batch_size: 15
interval_seconds: 45
objectives:
  - drug_discovery
  - green_chemistry
min_score_threshold: 0.4
max_memory_gb: 16.0
max_cpu_percent: 85.0
output_directory: "results/campaign_001"

seed_molecules:
  - "CCO"
  - "c1ccccc1" 
  - "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

target_properties:
  bioactivity: 0.8
  sustainability: 0.6
  reproducibility: 0.9
  efficiency: 0.7

validation_criteria:
  min_drug_likeness: 0.7
  min_safety_score: 0.8
  max_molecular_weight: 500
  max_lipinski_violations: 1
```

### AKG Configuration (YAML)

```yaml
# akg_config.yaml
neo4j:
  uri: "bolt://localhost:7687"
  user: "neo4j"
  password: "fotquantum"
  database: "neo4j"

schema:
  node_prefix: "FoTChem_"
  relationship_prefix: "FoTChem_"
  
export:
  format: "json"
  include_metadata: true
  max_discoveries: 1000
```

## 📊 Error Handling

### Exception Types

```python
class FoTChemistryError(Exception):
    """Base exception for FoTChemistry"""
    pass

class QuantumEngineError(FoTChemistryError):
    """Quantum engine initialization or operation failed"""
    pass

class MolecularGenerationError(FoTChemistryError):  
    """Molecular generation failed"""
    pass

class ValidationError(FoTChemistryError):
    """Molecular validation failed"""
    pass

class AKGConnectionError(FoTChemistryError):
    """AKG/Neo4j connection failed"""
    pass
```

### Error Handling Examples

```python
try:
    engine = ChemistryVQbitEngine(use_gpu=True)
except QuantumEngineError as e:
    logger.error(f"Quantum engine failed: {e}")
    # Fallback to CPU-only mode
    engine = ChemistryVQbitEngine(use_gpu=False)

try:
    candidates = generator.generate_candidates("CCO")
except MolecularGenerationError as e:
    logger.warning(f"Generation failed: {e}")
    candidates = []  # Continue with empty list

try:
    akg.store_discovery(discovery)
except AKGConnectionError as e:
    logger.error(f"Storage failed: {e}")
    # Fallback to local file storage
    with open("backup_discoveries.json", "a") as f:
        json.dump(discovery.__dict__, f)
```

## 🧪 Testing and Validation

### Unit Tests

```python
# test_quantum_engine.py
import unittest
from core.chemistry_vqbit_engine import ChemistryVQbitEngine, ChemistryPropertyType

class TestQuantumEngine(unittest.TestCase):
    def setUp(self):
        self.engine = ChemistryVQbitEngine(use_gpu=False)  # CPU for testing
    
    def test_vqbit_creation(self):
        scores = {ChemistryPropertyType.BIOACTIVITY: 0.5}
        vqbit = self.engine.create_molecular_vqbit(scores)
        self.assertEqual(len(vqbit.amplitudes), 8096)
        self.assertAlmostEqual(np.linalg.norm(vqbit.amplitudes), 1.0)
    
    def test_property_measurement(self):
        scores = {ChemistryPropertyType.BIOACTIVITY: 0.8}
        vqbit = self.engine.create_molecular_vqbit(scores)
        measurement = self.engine.measure_property(vqbit, ChemistryPropertyType.BIOACTIVITY)
        self.assertIsInstance(measurement, float)
        self.assertGreaterEqual(measurement, 0.0)
        self.assertLessEqual(measurement, 1.0)
```

### Integration Tests

```python
# test_discovery_pipeline.py
def test_end_to_end_discovery():
    """Test complete discovery pipeline"""
    config = ContinuousConfig(batch_size=3, objectives=["drug_discovery"])
    engine = ContinuousChemistryDiscoveryEngine(config)
    
    discoveries = engine.run_discovery_batch("test_batch")
    
    assert len(discoveries) <= 3
    for discovery in discoveries:
        assert discovery.combined_score >= config.min_score_threshold
        assert discovery.validation_passed
        assert discovery.smiles  # Valid SMILES
```

## 📚 Related Documentation

- **[Installation Guide](Installation-Guide)** - Setup and configuration
- **[Quick Start Tutorial](Quick-Start-Tutorial)** - Getting started
- **[Quantum vQbit Engine](Quantum-vQbit-Engine)** - Quantum mechanics details
- **[Molecular Discovery Pipeline](Molecular-Discovery-Pipeline)** - Discovery process
- **[Troubleshooting](Troubleshooting)** - Problem resolution

---

**🔧 Complete API reference for quantum-guided molecular discovery**
