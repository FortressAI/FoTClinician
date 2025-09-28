# Demo 3: Δ-Learning for Conformer Energies

## Overview

Demonstrate Δ-learning approach: baseline GFN2-xTB → Δ-learn to DFT on QM9-like subset. Generate truth claims about model accuracy and validate through reproducible pipelines with container images.

## Goals

- Baseline conformer energies with GFN2-xTB semiempirical method
- Δ-learn corrections to DFT level (B3LYP/6-31G*) on curated subset
- Generate claims: "Δ-model MAE ≤ X kcal/mol for molecule class Z"
- Independent reruns earn "replicated" badge for reproducibility

## Architecture

```
03_delta-learning-energies/
├── README.md                     # This file
├── run_demo.py                  # Main pipeline runner
├── data/
│   ├── qm9_subset.csv          # QM9 molecules subset
│   ├── conformers/             # Generated conformer geometries
│   └── energies/               # Computed energies
├── methods/
│   ├── xtb_baseline.py         # GFN2-xTB calculations
│   ├── dft_reference.py        # DFT reference calculations  
│   └── delta_learning.py       # Δ-learning ML model
├── pipelines/
│   ├── snakemake/             # Reproducible workflows
│   └── containers/            # Docker containers
└── claims/
    └── energy_accuracy_claims.jsonld
```

## Quick Start

```bash
# Generate conformers and compute energies
python methods/xtb_baseline.py --molecules data/qm9_subset.csv
python methods/dft_reference.py --subset data/qm9_subset.csv

# Train Δ-learning model
python methods/delta_learning.py --train

# Validate accuracy claims
python run_demo.py
```

## Expected Claims

- "GFN2-xTB + Δ-learning achieves MAE ≤ 2.5 kcal/mol on organic molecules"
- "Model generalizes to unseen scaffolds with MAE ≤ 4.0 kcal/mol"
- "Computational cost reduced 100× vs. direct DFT"

---

**Status**: 🏗️ Scaffold - Ready for Implementation  
**Contributors**: FoT Research Team
