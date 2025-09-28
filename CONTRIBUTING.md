# Contributing to FoTChemistry

Thank you for your interest in contributing to FoTChemistry! This project aims to create an open, truth-driven platform for collaborative chemistry research. We welcome contributions that maintain scientific integrity and advance the Field of Truth (FoT) methodology.

## 🎯 Quick Start

1. **Fork the repository** and clone your fork
2. **Set up the development environment** (see [Development Setup](#development-setup))
3. **Choose an issue** from our "good first issue" labels or propose a new contribution
4. **Follow our guidelines** for code quality, scientific integrity, and documentation
5. **Submit a pull request** with clear description and tests

## 🧪 Contributing Guidelines

### Scientific Integrity (Non-Negotiable)

All contributions must adhere to FoT principles:

#### ✅ **Required Practices**
- **No simulations or fake results** - All computational results must be real
- **Reproducible code** - Include environment specifications, seeds, versions
- **Evidence-based claims** - Support all assertions with appropriate data
- **Transparent limitations** - Document uncertainties and scope boundaries
- **Safety compliance** - Follow GHS guidelines and dual-use restrictions

#### ❌ **Rejected Contributions**
- Simulated or fabricated experimental data
- Exaggerated therapeutic claims without validation
- Hidden limitations or cherry-picked results
- Unsafe synthetic procedures without warnings
- Non-reproducible methods or environments

### Code Quality Standards

#### Python Code
```bash
# Format with black
black your_file.py

# Lint with flake8
flake8 your_file.py

# Type check with mypy
mypy your_file.py

# Test with pytest
pytest tests/
```

#### Chemistry-Specific Requirements
- **RDKit validation**: All molecules must pass `Chem.SanitizeMol()`
- **Unit consistency**: Use standard units (Kelvin, mol/L, kJ/mol)
- **Mass balance**: Chemical equations must balance within tolerance
- **Stereochemistry**: Properly handle chiral centers and isomers

#### Documentation Standards
- **Docstrings**: Use Google-style docstrings for all functions
- **Type hints**: Required for all function signatures
- **Examples**: Include usage examples in docstrings
- **Tests**: Unit tests required for all new functionality

## 🚀 Development Setup

### Prerequisites
```bash
# Python 3.9+ required
python --version

# Install dependencies
pip install -r requirements.txt

# Install development dependencies  
pip install -e ".[dev]"
```

### Docker Environment
```bash
# Start graph databases
cd akg && docker-compose up -d

# Verify services
docker-compose ps
```

### Pre-commit Hooks
```bash
# Install pre-commit
pip install pre-commit

# Set up hooks
pre-commit install

# Run on all files
pre-commit run --all-files
```

## 📁 Project Structure

Understanding our architecture helps you contribute effectively:

```
FoTChemistry/
├── ontology/           # OWL ontology and JSON-LD contexts
├── akg/               # Graph database infrastructure
├── agents/            # Truth-mining agents
│   ├── ingestion/     # Data ingestion from ChEMBL, PubChem, etc.
│   ├── validation/    # Chemical data validation
│   ├── measurement/   # Computational measurements
│   ├── ethics/        # Safety and ethics validation
│   └── collapse/      # Truth collapse algorithms
├── pipelines/         # Reproducible compute workflows
├── demos/             # Show-off projects and examples
└── datasets/          # Versioned datasets and examples
```

## 🎯 Contribution Types

### 1. **Good First Issues**

Perfect for new contributors:

- Map 50 RXNO terms into `fot_chemistry.owl`
- Write SPARQL query to flag mass-imbalanced reactions
- RDKit unit test: temperature/pressure parsing & normalization
- Add citation bot: ensure each claim has DOI/URL evidence fields
- Build tiny KG viewer (Streamlit/FastAPI) for claims & provenance

### 2. **Agent Development**

Contribute to our truth-mining agents:

```python
# Example: Validation Agent
class MassBalanceValidator(Agent):
    def validate_reaction(self, reaction: Reaction) -> ValidationResult:
        """Validate that reaction mass is balanced within tolerance."""
        # Your implementation here
        pass
```

### 3. **Ontology Enhancement**

Extend our chemistry ontology:

- Map to community ontologies (ChEBI, RXNO, CHMO)
- Add new chemical entity types
- Define vQbit quantum relationships
- Create SPARQL validation queries

### 4. **Demo Projects**

Contribute to our showcase demos:

- **Reaction-Yield Truth Mining**: ORD/USPTO data processing
- **pKa/logP Truth Claims**: Computational property prediction
- **Green-Chemistry Advisor**: Safer reaction alternatives
- **Delta-Learning Energies**: QM/MM energy predictions
- **Open Replication Challenge**: Community validation

### 5. **Pipeline Development**

Create reproducible computational workflows:

```yaml
# Example: Snakemake workflow for pKa prediction
rule predict_pka:
    input: "molecules.smi"
    output: "pka_predictions.csv"
    conda: "envs/rdkit.yaml"
    script: "scripts/predict_pka.py"
```

## 🔬 Testing Requirements

### Unit Tests
```python
def test_molecule_validation():
    """Test molecule sanitization and validation."""
    mol = Chem.MolFromSmiles("CCO")
    assert mol is not None
    assert Chem.SanitizeMol(mol) == Chem.SanitizeFlags.SANITIZE_NONE
```

### Integration Tests
```python
def test_truth_mining_workflow():
    """Test end-to-end truth mining workflow."""
    claim = create_test_claim()
    result = truth_mining_pipeline(claim)
    assert result.confidence > 0.5
    assert result.evidence is not None
```

### Property-Based Tests
```python
@given(smiles=valid_smiles_strategy())
def test_molecular_properties(smiles):
    """Test molecular property calculations are consistent."""
    mol = Chem.MolFromSmiles(smiles)
    mw1 = Descriptors.MolWt(mol)
    mw2 = rdMolDescriptors.CalcExactMolWt(mol)
    assert abs(mw1 - mw2) < 0.1  # Within tolerance
```

## 📊 Pull Request Process

### 1. **Preparation**
- Create feature branch: `git checkout -b feature/your-feature-name`
- Write tests first (TDD encouraged)
- Implement your changes
- Update documentation

### 2. **Quality Checks**
```bash
# Run all checks locally
make check-all

# Or individually:
black .
flake8 .
mypy .
pytest
```

### 3. **Pull Request Template**

Your PR should include:

```markdown
## Description
Brief description of changes and motivation.

## Type of Change
- [ ] Bug fix (non-breaking change)
- [ ] New feature (non-breaking change)
- [ ] Breaking change (fix/feature causing existing functionality to change)
- [ ] Documentation update
- [ ] Agent enhancement
- [ ] Ontology extension

## Scientific Validation
- [ ] All results are reproducible
- [ ] No simulated or fake data
- [ ] Limitations clearly documented
- [ ] Safety guidelines followed
- [ ] Evidence supports all claims

## Testing
- [ ] Unit tests pass
- [ ] Integration tests pass
- [ ] New tests added for new functionality
- [ ] Chemical validation tests included

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] Changes are backward compatible
```

### 4. **Review Process**

Your PR will be reviewed for:

- **Scientific accuracy**: Claims supported by evidence
- **Code quality**: Follows style guide and best practices
- **Safety compliance**: No hazardous procedures without warnings
- **Reproducibility**: Can others reproduce your results?
- **Documentation**: Clear and complete documentation

## 🌟 Recognition

Contributors are recognized through:

- **CRediT roles** in pull requests (conceptualization, methodology, software, etc.)
- **Contributor list** in README and documentation
- **ORCID integration** for academic recognition
- **Zenodo releases** with contributor metadata
- **Hall of Replications** for validation contributors

## 📞 Getting Help

### Community Channels
- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: General questions and brainstorming
- **Email**: research@fieldoftruth.org for private inquiries

### Weekly Community Calls
- **When**: Fridays at 3 PM UTC
- **Where**: Zoom link in GitHub discussions
- **What**: Demo progress, discuss issues, plan features

### Mentorship Program
New contributors can request mentorship for:
- Understanding FoT methodology
- Chemistry domain guidance
- Software architecture questions
- Publication pathways

## 🔒 Security and Safety

### Reporting Security Issues
- **Email**: security@fieldoftruth.org
- **Encryption**: PGP key available on request
- **Response time**: Within 72 hours

### Chemical Safety Guidelines
- No synthesis procedures for regulated compounds
- GHS hazard warnings required for dangerous chemicals
- Dual-use research guidelines strictly enforced
- Safety data sheets referenced when appropriate

## 📜 License and Copyright

By contributing, you agree that your contributions will be licensed under:
- **Code**: Apache License 2.0
- **Data**: CC-BY-4.0
- **Documentation**: CC-BY-4.0

### Contributor License Agreement
First-time contributors will be asked to sign a simple CLA confirming:
- You have the right to contribute the code
- You agree to the project licenses
- You understand the scientific integrity requirements

## 🎯 Roadmap Alignment

Current priorities (see full roadmap in README):

### **Week 1-2 (Foundation)**
- Ontology enhancements and mappings
- Graph database optimization
- CI/CD pipeline improvements

### **Week 3-4 (First Demo)**
- Reaction-yield truth mining demo
- Agent development and testing
- Documentation improvements

### **Week 5-6 (Launch)**
- Additional demos (pKa/logP, Green Chemistry)
- Zenodo integration
- Community building

Align your contributions with these priorities for maximum impact!

---

Thank you for contributing to FoTChemistry! Together, we're building the future of open, truth-driven chemistry research. 🧪✨
