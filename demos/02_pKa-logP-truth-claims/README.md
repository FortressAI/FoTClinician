# Demo 2: pKa & logP Truth Claims

## Overview

This demo demonstrates computational property prediction with truth-mining validation. We compute pKa and logP values for curated molecular sets using multiple computational methods (Psi4, OpenFF, RDKit), compare against experimental benchmarks, and establish truth claims through external replication.

## Goals

- Compute pKa/logP for curated molecular datasets
- Compare computational predictions against experimental benchmarks (SAMPL-like sets)
- Mark predictions as claims; collapse to truth when replicated within tolerance
- Demonstrate reproducible computational workflows with containerized environments

## Architecture

```
02_pKa-logP-truth-claims/
├── README.md                    # This file
├── run_demo.py                 # Main demo runner
├── data/                       # Input molecular datasets
│   ├── sampl_pka_set.csv      # SAMPL pKa challenge molecules
│   ├── drugbank_subset.csv    # DrugBank molecules with exp. data
│   ├── experimental/          # Experimental reference data
│   └── computed/              # Computational predictions
├── methods/                    # Computational methods
│   ├── rdkit_pka.py           # RDKit pKa prediction
│   ├── psi4_pka.py            # Psi4 quantum chemical pKa
│   ├── rdkit_logp.py          # RDKit logP prediction
│   └── openff_properties.py   # OpenFF force field properties
├── validation/                 # Truth validation agents
│   ├── experimental_validator.py  # Compare vs. experimental data
│   ├── method_validator.py        # Cross-method validation
│   ├── uncertainty_validator.py   # Uncertainty quantification
│   └── replication_validator.py   # External replication tracking
├── pipelines/                  # Computational workflows
│   ├── snakemake/             # Snakemake workflow definitions
│   ├── containers/            # Docker/Singularity containers
│   └── envs/                  # Conda environment specifications
├── claims/                     # Generated truth claims
│   ├── pka_claims.jsonld      # pKa truth claims in JSON-LD
│   ├── logp_claims.jsonld     # logP truth claims in JSON-LD
│   └── collapsed_truth.jsonld # Validated truth claims
├── notebooks/                  # Analysis notebooks
│   ├── 01_method_comparison.ipynb
│   ├── 02_uncertainty_analysis.ipynb
│   └── 03_truth_validation.ipynb
└── tests/                      # Unit tests
    ├── test_pka_methods.py
    ├── test_logp_methods.py
    └── test_validation.py
```

## Quick Start

```bash
# Install computational chemistry dependencies
conda env create -f pipelines/envs/comp_chem.yaml
conda activate comp_chem

# Run pKa predictions
python methods/rdkit_pka.py --input data/sampl_pka_set.csv
python methods/psi4_pka.py --input data/sampl_pka_set.csv

# Run validation and truth mining
python run_demo.py

# Explore results
jupyter lab notebooks/01_method_comparison.ipynb
```

## Computational Methods

### pKa Prediction Methods

#### RDKit pKa (Fast, Empirical)
```python
from rdkit import Chem
from rdkit.Chem import Descriptors

def predict_pka_rdkit(smiles: str) -> Dict[str, float]:
    """Predict pKa using RDKit empirical methods."""
    mol = Chem.MolFromSmiles(smiles)
    # Implementation using RDKit descriptors
    return {"pka_acid": pka_value, "uncertainty": std_dev}
```

#### Psi4 pKa (Accurate, Quantum Chemical)
```python
import psi4

def predict_pka_psi4(smiles: str, method: str = "B3LYP/6-31G*") -> Dict[str, float]:
    """Predict pKa using DFT calculations in Psi4."""
    # Generate 3D conformer
    # Optimize geometry 
    # Calculate solvation free energies
    # Compute thermodynamic cycle
    return {"pka_acid": pka_value, "uncertainty": std_dev}
```

### logP Prediction Methods

#### RDKit logP (Multiple Descriptors)
```python
def predict_logp_rdkit(smiles: str) -> Dict[str, float]:
    """Predict logP using multiple RDKit methods."""
    mol = Chem.MolFromSmiles(smiles)
    methods = {
        "crippen": Crippen.MolLogP(mol),
        "wildman": rdMolDescriptors.CalcWildmanKrilovLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol)
    }
    return {"logp_consensus": consensus_value, "uncertainty": std_dev}
```

## Truth Claims Framework

### Claim Types

#### Individual Method Claims
```json
{
  "@context": "../../ontology/fot_chemistry_context.jsonld",
  "@type": "Claim", 
  "@id": "claim:pka-rdkit-aspirin",
  "hasConfidence": 0.75,
  "hasVirtueWeight": 0.68,
  "label": "RDKit predicts aspirin pKa = 3.48 ± 0.3",
  "hasMeasurement": {
    "@type": "Measurement",
    "method": "RDKit empirical pKa",
    "value": 3.48,
    "uncertainty": 0.3,
    "unit": "pH_units"
  },
  "hasEvidence": {
    "@type": "Evidence", 
    "experimental_reference": 3.5,
    "literature_doi": "10.1021/acs.jcim.example"
  }
}
```

#### Consensus Claims (Higher Confidence)
```json
{
  "@context": "../../ontology/fot_chemistry_context.jsonld",
  "@type": "Claim",
  "@id": "claim:pka-consensus-aspirin", 
  "hasConfidence": 0.91,
  "hasVirtueWeight": 0.87,
  "label": "Consensus pKa for aspirin = 3.52 ± 0.15",
  "hasMeasurement": {
    "@type": "Measurement",
    "method": "Multi-method consensus",
    "value": 3.52,
    "uncertainty": 0.15,
    "contributing_methods": ["RDKit", "Psi4", "OpenFF"]
  }
}
```

### Truth Collapse Criteria

Claims are collapsed to accepted truth when:

1. **Experimental Agreement**: Prediction within 0.5 pKa units or 0.3 logP units of experiment
2. **Method Consensus**: At least 3 different methods agree within uncertainty
3. **External Replication**: Independent group reproduces prediction within tolerance  
4. **Uncertainty Quantified**: Proper error bars and confidence intervals provided

## Validation Agents

### Experimental Validator
- Compare predictions against SAMPL challenge experimental data
- Flag predictions outside acceptable error bounds
- Track method-specific systematic biases

### Method Validator  
- Cross-validate predictions between different computational methods
- Identify outliers and method-specific failures
- Compute consensus values with uncertainty propagation

### Uncertainty Validator
- Ensure all predictions include proper uncertainty estimates
- Validate that confidence intervals contain experimental values
- Flag overconfident predictions

### Replication Validator
- Track external replications of computational predictions
- Monitor for computational environment differences
- Validate reproducibility across different hardware/software

## Expected Results

### Success Metrics
- **Coverage**: 100+ molecules with pKa and logP predictions
- **Accuracy**: 80% of predictions within error bounds of experiment
- **Consensus**: 60% of molecules have method consensus within 0.3 units
- **Replications**: 10+ external replications confirming predictions

### Benchmark Comparisons
- **SAMPL pKa**: Compare against SAMPL challenge experimental data
- **DrugBank**: Validate against pharmaceutical compound properties  
- **Literature**: Cross-reference with published experimental values

## Reproducible Workflows

### Snakemake Pipeline
```python
# Snakefile for pKa/logP prediction pipeline
rule all:
    input: "claims/collapsed_truth.jsonld"

rule rdkit_pka:
    input: "data/{dataset}.csv"
    output: "computed/rdkit_pka_{dataset}.csv"
    conda: "envs/rdkit.yaml"
    script: "methods/rdkit_pka.py"

rule psi4_pka:
    input: "data/{dataset}.csv"  
    output: "computed/psi4_pka_{dataset}.csv"
    conda: "envs/psi4.yaml"
    resources: cpu=4, mem_mb=8000
    script: "methods/psi4_pka.py"

rule validate_claims:
    input: 
        rdkit="computed/rdkit_pka_{dataset}.csv",
        psi4="computed/psi4_pka_{dataset}.csv"
    output: "claims/pka_claims_{dataset}.jsonld"
    script: "validation/generate_claims.py"
```

### Container Definitions
```dockerfile
# Dockerfile for Psi4 quantum chemistry
FROM psi4/psi4:latest

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY methods/ /app/methods/
WORKDIR /app

ENTRYPOINT ["python", "methods/psi4_pka.py"]
```

## Contributing

### Enhancement Opportunities
- Add new computational methods (COSMO-RS, ChemAxon)
- Implement machine learning property prediction models
- Expand to additional molecular properties (solubility, permeability)
- Create web interface for property prediction and validation

### Good First Issues
- [ ] Implement SMILES validation and canonicalization
- [ ] Add unit tests for uncertainty propagation
- [ ] Create visualization of method agreement/disagreement
- [ ] Implement automated literature search for experimental values
- [ ] Add support for ionic species and protonation states

## References

- **SAMPL Challenges**: https://samplchallenges.github.io/
- **RDKit Documentation**: https://www.rdkit.org/docs/
- **Psi4 Manual**: https://psicode.org/psi4manual/master/
- **OpenFF Toolkit**: https://open-forcefield-toolkit.readthedocs.io/

---

**Status**: 🚧 Under Development  
**Last Updated**: September 2025  
**Contributors**: FoT Research Team
