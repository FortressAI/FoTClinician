# Demo 1: Reaction-Yield Truth Mining (USPTO → ORD)

## Overview

This demo showcases truth-mining workflows on reaction yield data from the USPTO and Open Reaction Database (ORD). Agents propose claims about catalyst effects on reaction yields, which are then validated through logic checks, statistical analysis, and reproducibility assessment.

## Goals

- Ingest open reaction records from ORD/USPTO sample
- Generate testable claims: "Under conditions X, catalyst Y improves yield Δ%"
- Validate claims through multiple validation agents
- Output signed claims with evidence bundles for Zenodo deposit

## Architecture

```
01_reaction-yield-truth-mining/
├── README.md                 # This file
├── run_demo.py              # Main demo runner
├── data/                    # Input datasets
│   ├── uspto_sample.json    # Sample USPTO reaction data
│   ├── ord_sample.json      # Sample ORD reaction data
│   └── processed/           # Processed datasets
├── agents/                  # Demo-specific agents
│   ├── ingestion_agent.py   # Data ingestion from ORD/USPTO
│   ├── claim_generator.py   # Generate yield improvement claims
│   ├── logic_validator.py   # Mass balance, unit checks
│   ├── stats_validator.py   # Statistical significance testing
│   └── repro_validator.py   # Reproducibility assessment
├── pipelines/               # Processing workflows
│   ├── snakemake/          # Snakemake workflows
│   └── cwl/                # CWL workflows  
├── outputs/                 # Generated claims and evidence
│   ├── claims/             # Validated claims in JSON-LD
│   ├── evidence/           # Supporting evidence bundles
│   └── zenodo/             # Zenodo export packages
├── notebooks/              # Jupyter analysis notebooks
│   ├── 01_data_exploration.ipynb
│   ├── 02_claim_generation.ipynb
│   └── 03_validation_analysis.ipynb
└── tests/                  # Unit tests for demo components
    ├── test_ingestion.py
    ├── test_validation.py
    └── test_claims.py
```

## Quick Start

```bash
# Install dependencies
pip install -r ../../requirements.txt

# Run the demo
python run_demo.py

# Explore in Jupyter
jupyter lab notebooks/01_data_exploration.ipynb
```

## Validation Criteria

### Logic Checks (Logic Validator)
- **Mass balance**: Reactants → Products must conserve mass within 1%
- **Units consistency**: Temperature (K), pressure (bar), concentration (mol/L)
- **Temperature/pressure plausibility**: Realistic experimental conditions
- **Catalyst stoichiometry**: Reasonable catalyst loading (0.1-20 mol%)

### Statistical Checks (Stats Validator)
- **Effect size**: Yield improvement > 5% absolute difference
- **Confidence intervals**: 95% CI for yield difference
- **P-values**: Significance test p < 0.05 for improvement claims
- **Bayesian evidence ratio**: Bayes factor > 3 for strong evidence

### Reproducibility Flags (Repro Validator)
- **Dataset diversity**: Multiple reaction classes represented
- **Replicates present**: At least 3 independent measurements
- **Experimental details**: Complete reaction conditions reported
- **Control experiments**: Appropriate baseline reactions included

## Sample Claims Generated

### Example Claim 1: Catalyst Effect
```json
{
  "@context": "../../ontology/fot_chemistry_context.jsonld",
  "@type": "Claim",
  "@id": "claim:catalyst-effect-001",
  "hasConfidence": 0.87,
  "hasVirtueWeight": 0.92,
  "label": "Pd(PPh3)4 improves Suzuki coupling yield by 15.3%",
  "hasEvidence": [
    {
      "@type": "Evidence",
      "description": "Statistical analysis of 47 Suzuki reactions",
      "confidence_interval": [12.1, 18.5],
      "p_value": 0.003,
      "bayes_factor": 8.4
    }
  ],
  "validatedBy": ["LogicValidator", "StatsValidator", "ReproValidator"]
}
```

### Example Claim 2: Solvent Effect
```json
{
  "@context": "../../ontology/fot_chemistry_context.jsonld", 
  "@type": "Claim",
  "@id": "claim:solvent-effect-001",
  "hasConfidence": 0.73,
  "hasVirtueWeight": 0.81,
  "label": "DMF solvent increases amide formation yield by 22%",
  "hasEvidence": [
    {
      "@type": "Evidence",
      "description": "Meta-analysis of 83 amide coupling reactions",
      "confidence_interval": [18.2, 25.8],
      "p_value": 0.001
    }
  ]
}
```

## Data Sources

### USPTO Dataset Sample
- **Source**: USPTO reaction data via ChemRxiv
- **Size**: 1,000 representative reactions
- **Format**: JSON with SMILES, conditions, yields
- **License**: Public domain

### ORD Dataset Sample  
- **Source**: Open Reaction Database
- **Size**: 500 well-documented reactions
- **Format**: Protocol buffer → JSON conversion
- **License**: CC-BY-4.0

## Truth Mining Workflow

1. **Ingestion**: Parse USPTO/ORD data into standardized format
2. **Claim Generation**: Identify potential catalyst/solvent effects
3. **Logic Validation**: Check mass balance and unit consistency
4. **Statistical Validation**: Test significance of claimed effects
5. **Reproducibility Assessment**: Evaluate experimental rigor
6. **Truth Collapse**: Collapse validated claims to accepted truth
7. **Evidence Packaging**: Bundle claims with supporting evidence
8. **Zenodo Export**: One-click export for permanent archival

## Expected Outcomes

### Success Metrics
- **Claims Generated**: 50-100 testable claims about reaction conditions
- **Validation Rate**: 60-80% of claims pass all validation criteria  
- **Truth Collapsed**: 20-30 high-confidence claims accepted as truth
- **Zenodo DOI**: Permanent archive with searchable metadata

### Failure Modes
- **Poor data quality**: Missing reaction conditions or yields
- **Statistical insignificance**: Effects too small to detect reliably
- **Reproducibility issues**: Insufficient experimental details
- **Mass balance violations**: Inconsistent stoichiometry

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for general guidelines.

### Demo-Specific Contributions
- Add new validation criteria for chemical feasibility
- Implement machine learning models for yield prediction
- Expand to additional reaction databases (Reaxys, SciFinder)
- Create interactive visualization of validation results

### Good First Issues
- [ ] Implement temperature range validation (273-473 K reasonable)
- [ ] Add molecular weight calculation for mass balance
- [ ] Create unit tests for SMILES parsing and validation
- [ ] Build Streamlit dashboard for claim exploration
- [ ] Add export functionality for different file formats

---

**Status**: 🚧 Under Development  
**Last Updated**: September 2025  
**Contributors**: FoT Research Team
