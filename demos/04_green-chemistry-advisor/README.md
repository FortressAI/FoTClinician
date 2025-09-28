# Demo 4: Green-Chemistry Advisor

## Overview

Given a reaction plan, agents suggest safer solvent and better E-factor alternatives from open data (CHEM21/Glaxo guidelines). Outputs are advice claims (not procedures), gated by ethics agents for safety.

## Goals

- Parse reaction plans and identify improvement opportunities
- Suggest safer solvents based on CHEM21 solvent selection guide
- Recommend process improvements for better atom economy and E-factor
- Generate safety-validated advice claims with proper disclaimers

## Architecture

```
04_green-chemistry-advisor/
├── README.md                    # This file
├── run_demo.py                 # Main advisor system
├── data/
│   ├── chem21_solvents.csv     # CHEM21 solvent safety data
│   ├── reaction_examples.json  # Example reaction plans
│   └── safety_guidelines/      # GHS and safety data
├── advisors/
│   ├── solvent_advisor.py      # Solvent replacement suggestions
│   ├── process_advisor.py      # Process improvement suggestions
│   └── safety_validator.py     # Safety and ethics validation
├── knowledge/
│   ├── solvent_db.py          # Solvent properties database
│   ├── safety_rules.py        # Safety rule engine
│   └── green_metrics.py       # Green chemistry metrics
└── advice_claims/
    └── green_advice_claims.jsonld
```

## Quick Start

```bash
# Load reaction plan
python run_demo.py --reaction "data/reaction_examples.json"

# Get green chemistry advice  
python advisors/solvent_advisor.py --input "CCCl3" --output suggestions.json

# Validate safety compliance
python advisors/safety_validator.py --check suggestions.json
```

## Example Advice Claims

- "Replace dichloromethane with ethyl acetate for improved safety profile"
- "Consider microwave heating to reduce reaction time by 75%"
- "Atom economy can be improved from 45% to 78% by eliminating protecting group"

## Safety Guardrails

- No synthesis procedures for regulated/hazardous compounds
- All suggestions include proper safety warnings
- Ethics agent validates all advice claims before publication

---

**Status**: 🏗️ Scaffold - Ready for Implementation  
**Contributors**: FoT Research Team
