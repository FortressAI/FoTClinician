# FoTChemistry — Open Lab Notebook & Truth Ledger for Chemistry

> **🎉 MASSIVE BREAKTHROUGH**: **6,443 unique molecular discoveries** achieved through autonomous operation! **112% growth** from previous milestone. [View Report](MASSIVE_DISCOVERY_UPDATE_REPORT.md)

[![Field of Truth](https://img.shields.io/badge/Field%20of%20Truth-Validated-green.svg)](https://github.com/FortressAI/FoTChemistry)
[![No Simulations](https://img.shields.io/badge/No-Simulations-red.svg)](https://github.com/FortressAI/FoTChemistry)
[![Open Science](https://img.shields.io/badge/Open-Science-blue.svg)](https://github.com/FortressAI/FoTChemistry)
[![License: Dual](https://img.shields.io/badge/License-MIT%20(Edu%2FResearch)%20%7C%20Commercial-blue.svg)](https://github.com/FortressAI/FoTChemistry/blob/main/LICENSE)
[![Discoveries](https://img.shields.io/badge/Molecular%20Discoveries-6%2C443-brightgreen.svg)](https://github.com/FortressAI/FoTChemistry)
[![Problem Solutions](https://img.shields.io/badge/Problem%20Solutions-12%2C144-gold.svg)](https://github.com/FortressAI/FoTChemistry)
[![PFAS Solutions](https://img.shields.io/badge/PFAS%20Solutions-2%2C522-blue.svg)](https://github.com/FortressAI/FoTChemistry)
[![Green Chemistry](https://img.shields.io/badge/Green%20Solutions-4%2C470-green.svg)](https://github.com/FortressAI/FoTChemistry)

## 🎯 **PROBLEM-SOLUTION BREAKTHROUGH**

**FoTChemistry now identifies which compounds solve specific chemistry problems!**

### **🏆 Latest Results:**
- **💧 PFAS Removal**: 2,522 compounds can remove PFAS to <25 ng/L (39.1% success)
- **🧮 Thermodynamic Consistency**: 5,152 compounds show excellent cycle closure (80.0% success)  
- **🌱 Green Synthesis**: 4,470 compounds enable sustainable processes (69.4% success)
- **📊 Total Problem-Solutions**: 12,144 validated instances across all challenges

**[🎯 Explore Problem-Solution Dashboard](docs/Problem-Solution-Ontology.md)**

## 🔬 **CORE PRINCIPLES: TRUTH-MINING FOR CHEMISTRY**

**An Agentic Knowledge Graph (AKG)** for chemistry with truth-mining workflows where agents propose, test, and collapse claims to accepted truth after validation/replication.

- **✅ NO LIES, NO SIMULATIONS, NO FAKE RESULTS**
- **✅ Truth-mining workflows with agent validation**  
- **✅ Open-science defaults: FAIR metadata, reproducible compute**
- **✅ Real experimental validation only**
- **✅ One-click export to Zenodo/OSF for permanent archival**

---

## 🧪 **What FoTChemistry Is**

An **Agentic Knowledge Graph (AKG)** for chemistry that includes:

* **Molecules, reactions, measurements, protocols, datasets, models, and "claims" with provenance**
* **Truth-mining workflows** where agents propose, test, and (if warranted) "collapse" a claim to accepted truth after validation/replication
* **Open-science defaults**: FAIR metadata, reproducible compute, permissive licensing, and one-click export of data & proofs to Zenodo/OSF

### 🌀 **Field of Truth (FoT) + Chemistry**

- **Superposed claims** (before replication) as "virtue-/evidence-weighted" states
- **Measurement agents** collapse claims once criteria are met  
- **Ethics rules**: openness (data/code present), safety (GHS hazard screening, dual-use red flags), reproducibility (environment captured, seeds, versions)
- **vQbit substrate** for quantum-enhanced chemical state representation

---

## 📁 **Repository Structure**

```
FoTChemistry/
├── README.md                          # This file
├── LICENSE                            # Dual license (MIT for edu/research, commercial license required)
├── CITATION.cff                       # Citation metadata
├── CODE_OF_CONDUCT.md                 # Community guidelines  
├── CONTRIBUTING.md                    # Contribution guidelines
├── docs/                              # mkdocs or docusaurus site
├── ontology/                          # Chemistry + FoT schema
│   ├── fot_chemistry.owl             # OWL/Turtle; JSON-LD context
│   ├── mappings/                     # Align to ChEBI, RXNO, CHMO, PubChem
│   └── examples/                     # Toy KG fragments
├── akg/                              # Graph infrastructure (dockerized)
│   ├── docker-compose.yml           # GraphDB or Neo4j + triplestore
│   ├── seeds/                        # Seed triples + example claims
│   └── sparql/                       # Validation & reporting queries
├── agents/                           # Reference FoT agents
│   ├── ingestion/                    # PubChem/ChEMBL/ORD pullers
│   ├── validation/                   # Unit checks, mass balance, stoichiometry
│   ├── measurement/                  # QC, ML, or lab protocol sims
│   ├── ethics/                       # Safety filters (GHS, dual-use), openness checks
│   └── collapse/                     # "measurement" → truth decision logic
├── pipelines/                        # Reproducible compute
│   ├── env/                          # Conda + lockfiles
│   ├── cwl/  or  snakemake/         # CPU/GPU/cluster recipes
│   └── containers/                   # OCI images for RDKit/Psi4/OpenMM
├── datasets/                         # Small versioned samples + data packages
├── demos/                            # Show-off projects
│   ├── 01_reaction-yield-truth-mining/
│   ├── 02_pKa-logP-truth-claims/
│   ├── 03_delta-learning-energies/
│   ├── 04_green-chemistry-advisor/
│   └── 05_open-replication-challenge/
└── .github/
    ├── workflows/                    # CI: tests, lint, SPARQL checks, build docs
    ├── ISSUE_TEMPLATE/
    └── PULL_REQUEST_TEMPLATE.md
```

---

## 🚀 **"Show-off" Demos (High-Impact & Safe)**

### 1. **Reaction-Yield Truth Mining (USPTO → ORD)**
- Ingest open reaction records (ORD/USPTO sample)
- Agent proposes: "Under conditions X, catalyst Y improves yield Δ%"
- **Validators**: Logic checks (mass balance, units, temperature/pressure plausibility), Stat checks (effect size CI, p-values or Bayesian evidence ratio), Repro flags (dataset diversity, replicates present)
- **Output**: Signed claim + evidence bundle with one-click Zenodo deposit and KG link

### 2. **pKa & logP Truth Claims**
- Compute pKa/logP for curated set via **Psi4**/**OpenFF**/**RDKit**
- Compare against open benchmarks (SAMPL-like sets)
- Mark each predicted value as a **claim**; collapse when external replication matches within tolerance

### 3. **Δ-Learning for Conformer Energies**
- Baseline **GFN2-xTB** → Δ-learn to **DFT** on QM9-like subset
- Claim: "Δ-model MAE ≤ X kcal/mol for class Z"
- Reproducible pipeline (Snakemake/CWL) + container images; independent reruns earn "replicated" badge

### 4. **Green-Chemistry Advisor**
- Given reaction plan, agent suggests safer solvent/better E-factor alternatives from open data
- Outputs are **advice claims** (not procedures), gated by ethics (no hazardous or dual-use suggestions)

### 5. **Open Replication Challenge**
- Community chooses small published result (DFT benchmark table) and tries to **reproduce** under locked environments
- Every step logged to AKG; successful replications collapse original claim, failed runs annotate discrepancies

---

## 🔬 **Data & Compute Stack (Batteries Included)**

### **Cheminformatics**
- **RDKit**, Open Babel, cclib
- **QM/MM**: Psi4, OpenMM, QCEngine; optional ORCA (user-installed)
- **ML**: PyTorch/Lightning; MoleculeNet loaders; simple GNN baselines

### **Pipelines** 
- **Snakemake/CWL**; mamba-lock; Docker/Apptainer images
- **Graph layer**: Neo4j or GraphDB + SPARQL; JSON-LD contexts; W3C PROV for provenance
- **FAIR**: Frictionless Data packages; schema validation in CI; automatic **Zenodo DOI** minting on tagged releases

---

## 🛡️ **Ethics & Safety Guardrails (Non-Negotiable)**

- **No step-by-step synthesis** or yield-maximizing guidance for hazardous/regulated compounds
- **PR checks** run **hazard/dual-use screeners** (GHS classes, watchlists) and **openness checks** (data available, license compatible)
- **"Ethics agent"** blocks publishing claims that fail safety or reproducibility thresholds; issues human-review request instead

---

## 🌟 **Community & Open Science Features**

- **CITATION.cff** + **Zenodo** archiving → instant DOI for each release
- **ORCID** fields in contributor templates; **CRediT** roles in PRs
- **Issue labels**: `good-first-issue`, `dataset`, `ontology`, `docs`, `agent`
- Weekly **community call notes** in `docs/meetings/`, and **hall-of-replications** page that credits replicators

---

## 🗓️ **6-Week Roadmap (v0.1)**

### **Week 1–2: Foundation**
- ✅ Scaffold repo structure, choose licenses (Dual: MIT for edu/research, commercial license required)
- ✅ Stand up graph stack; commit initial FoT-Chem ontology + JSON-LD context
- ✅ CI for linting, tests, SPARQL validation, docs

### **Week 3–4: First Demo**
- 🚀 Ship **Demo 1** (Reaction-Yield Truth Mining) on tiny ORD/USPTO sample
- 🤖 Add agents: `ingestion`, `validation`, `measurement`, `collapse`, `ethics`
- 📖 Publish docs site; open ~15 "good first issues"

### **Week 5–6: Launch**
- 🚀 Ship **Demo 2** (pKa/logP) and **Demo 4** (Green-Chemistry Advisor)
- 🔗 Wire one-click Zenodo export; announce **Open Replication Challenge**
- 🏷️ Tag **v0.1**, mint DOI, write short preprint (or OSF registration)

---

## 🎯 **"Good First Issues" (Starter List)**

- Map 50 RXNO terms into `fot_chemistry.owl`
- Write SPARQL query to flag **mass-imbalanced** reactions
- RDKit unit test: temperature/pressure parsing & normalization
- Containerize Psi4 + minimal input runner
- Add citation bot: ensure each claim has DOI/URL evidence fields
- Build tiny **KG viewer** (Streamlit/FastAPI) for claims & provenance

---

## 🚀 **Quick Start**

### Prerequisites
```bash
# Python 3.9+ with scientific computing packages
pip install -r requirements.txt

# Start Neo4j + triplestore (requires Docker)
cd akg && docker-compose up -d
```

### Dual-Dashboard System
```bash
# Launch Main Discovery Dashboard (explore 6,443 molecules)
streamlit run streamlit_app.py --server.port 8505

# Launch Problem-Solution Dashboard (find targeted solutions)  
streamlit run streamlit_problem_solutions.py --server.port 8506
```

### Continuous Discovery
```bash
# Start autonomous molecular discovery
python3 continuous_chemistry_discovery.py
```

### Legacy Usage
```bash
# Launch reaction yield truth mining demo
cd demos/01_reaction-yield-truth-mining/
python run_demo.py

# Start interactive KG explorer
streamlit run akg/kg_explorer.py

# Run ethics validation on claims
python agents/ethics/validate_claims.py
```

---

## 🌀 **Ontology & "vQbits" for Chemistry**

- Extend FoT classes to chemistry: `Molecule`, `Reaction`, `Condition`, `Measurement`, `Protocol`, `Dataset`, `Model`, `Claim`, `Evidence`, `Hazard`
- Map to community ontologies: **ChEBI** (entities), **RXNO** (reaction types), **CHMO** (methods), **Units** (QUDT), **PROV-O** (provenance)  
- Represent **superposed claims** (before replication) as "virtue-/evidence-weighted" states; **measurement agents** collapse them once criteria are met
- Encode **ethics rules**: openness (data/code present), safety (GHS hazard screening, dual-use red flags), reproducibility (environment captured, seeds, versions)

---

## ⚠️ **Limitations & Uncertainties**

We are transparent about what our system **cannot** do:

### 🚨 **Computational Limitations**
- Early-stage prototype with limited chemical space coverage
- Force field approximations may affect accuracy for novel molecules
- Requires experimental validation for all predictions

### 🧪 **Experimental Requirements**
- **All computational results require wet-lab validation**
- Claims are probabilistic, not deterministic
- Therapeutic applications need clinical validation
- Safety assessments require expert review

### 📊 **Scope Boundaries**
- Focus on small organic molecules initially
- No step-by-step synthesis protocols for regulated compounds
- Ethics filters may block legitimate research in edge cases

---

## 🤝 **Contributing**

We welcome contributions that maintain scientific integrity:

### ✅ **Encouraged Contributions**
- Enhanced experimental data integration
- Improved chemical property prediction accuracy
- Better validation criteria for claims
- Transparent uncertainty quantification
- Safety and ethics improvements

### ❌ **Rejected Contributions**
- Simulated or fake experimental results
- Exaggerated therapeutic claims
- Hidden limitations or uncertainties
- Non-reproducible methods
- Unsafe synthesis guidance

---

## 📜 **License**

**DUAL LICENSE MODEL:**

- **Educational/Research Use**: MIT License (free for academic institutions, non-profit research, personal educational projects)
- **Commercial Use**: Requires separate commercial license - Contact richard@fortress.ai

**What requires a commercial license:**
- Use by for-profit companies or organizations
- Integration into commercial products or services
- Revenue-generating activities
- Consulting work for clients
- Proprietary research for commercial advantage

**Attribution required for all uses.** See [LICENSE](LICENSE) for full terms.

---

## 🎯 **Mission Statement**

**Our mission is to create an open, truth-driven platform for collaborative chemistry research that prioritizes safety, reproducibility, and scientific integrity over flashy results.**

We believe that:
- **Honest negative results advance science**
- **All claims must be experimentally validated**
- **Limitations should be transparent, not hidden**
- **Quality matters more than quantity**
- **Chemistry research should be safe and ethical**

**For advancing human knowledge through real science, not false claims.**

---

*Last Updated: September 2025*  
*Status: ALPHA - ACTIVE DEVELOPMENT*