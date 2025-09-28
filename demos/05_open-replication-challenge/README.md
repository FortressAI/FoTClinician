# Demo 5: Open Replication Challenge

## Overview

Community chooses a small published result (e.g., DFT benchmark table) and attempts to reproduce it under locked computational environments. Every step is logged to the AKG; successful replications collapse the original claim, failed runs annotate discrepancies.

## Goals

- Select community-nominated published computational results for replication
- Provide locked computational environments (containers) for reproducibility
- Track all replication attempts in knowledge graph with full provenance
- Award "replicated" badges for successful reproductions

## Architecture

```
05_open-replication-challenge/
├── README.md                    # This file
├── challenges/
│   ├── challenge_001/          # First replication challenge
│   │   ├── original_paper.pdf  # Reference publication
│   │   ├── environment/        # Locked computational environment
│   │   ├── replication_attempts/ # Community attempts
│   │   └── results.jsonld      # Aggregated results
│   └── challenge_template/     # Template for new challenges
├── infrastructure/
│   ├── challenge_manager.py    # Challenge lifecycle management
│   ├── environment_builder.py  # Container environment setup
│   └── result_aggregator.py    # Collect and validate results
├── validation/
│   ├── numerical_validator.py  # Numerical agreement validation
│   ├── method_validator.py     # Methodology compliance
│   └── provenance_tracker.py   # Full provenance logging
└── community/
    ├── leaderboard.py          # Community contribution tracking
    └── badges.py               # Replication badge system
```

## Quick Start

```bash
# Browse active challenges
python challenges/list_active.py

# Join a challenge
python infrastructure/challenge_manager.py --join challenge_001

# Run replication attempt
docker run --rm -v $(pwd):/workspace fotchem/challenge-001:latest

# Submit results  
python infrastructure/result_aggregator.py --submit results.json
```

## Challenge Types

### Computational Benchmarks
- DFT energy calculations on small molecule sets
- Molecular property predictions (logP, pKa, solubility)
- Reaction barrier calculations
- Conformational analysis results

### Replication Criteria
- **Numerical Agreement**: Within 5% or 2σ of original results
- **Method Compliance**: Same computational methods and parameters
- **Environment Reproducibility**: Results match across different hardware
- **Full Provenance**: Complete audit trail of computational steps

## Community Features

- **Leaderboard**: Track successful replications by contributor
- **Badges**: "Replicator", "Method Expert", "Environment Builder"
- **Hall of Replications**: Showcase successful community validations
- **Collaboration**: Comment and discuss replication challenges

---

**Status**: 🏗️ Scaffold - Ready for Implementation  
**Contributors**: FoT Research Team
