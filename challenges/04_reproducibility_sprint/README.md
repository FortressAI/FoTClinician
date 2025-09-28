# Reproducibility Sprint: Truth-Mining Literature Claims

**Challenge Period**: Monthly events (ongoing)  
**Scientific Impact**: Restore trust in chemical literature  
**Target**: `≥20 claims collapsed/quarter; full provenance`

## 🎯 Sprint Mission

Transform disputed literature results into validated "claims with evidence" through crowd-replication under locked computational environments and collapse to truth using FoT validation rules.

## 📊 Success Metrics

| Metric | Monthly Target | Quarterly Target | Annual Target |
|--------|---------------|------------------|---------------|
| **Literature Claims Tested** | ≥7 | ≥20 | ≥80 |
| **Independent Replications** | ≥3 per claim | ≥60 total | ≥240 total |
| **Truth Collapses** | ≥5 | ≥15 | ≥60 |
| **Disputed Claims Resolved** | ≥2 | ≥6 | ≥25 |
| **Community Participants** | ≥10 | ≥30 | ≥100 |

## 🎯 Monthly Sprint Format

### **Sprint Structure (4 Weeks)**
```yaml
week_1_nomination:
  - community_nominates_claims
  - expert_panel_selects_targets
  - preparation_of_test_environments
  - protocol_standardization

week_2_3_replication:
  - distributed_replication_attempts
  - real_time_progress_tracking
  - discussion_and_troubleshooting
  - preliminary_results_sharing

week_4_validation:
  - results_compilation_and_analysis
  - statistical_evaluation
  - truth_collapse_decisions
  - recognition_and_documentation
```

### **Claim Selection Criteria**
```yaml
priority_targets:
  impact: "high_citation_count OR controversial"
  reproducibility: "conflicting_reports OR failed_replications"
  scope: "computationally_feasible"
  safety: "no_hazardous_procedures"
  
example_categories:
  - pka_logp_predictions
  - conformer_energies
  - reaction_barriers
  - catalyst_activities
  - solubility_predictions
  - binding_affinities
```

## 🔬 Current Sprint Targets

### **Sprint #1 (October 2025): pKa Prediction Methods**
```yaml
target_claim: "ChemAxon pKa predictor achieves MAE <0.3 pH units on SAMPL dataset"
original_source: "doi:10.1021/acs.jcim.example"
controversy: "multiple_groups_report_conflicting_accuracies"
test_scope: "100_molecule_subset_from_SAMPL6"
expected_participants: 15
computational_requirements: "minimal_ChemAxon_license_or_alternatives"
```

#### **Validation Protocol**
```yaml
dataset:
  molecules: "SAMPL6_pKa_challenge_subset_100"
  experimental_values: "literature_curated_high_confidence"
  molecule_diversity: "drug_like_organic_acids_bases"

computational_setup:
  methods_to_test: ["ChemAxon", "RDKit", "Epik", "ACD_Labs"]
  containers: "docker_images_with_locked_versions"
  environment: "ubuntu_20.04_python_3.9"
  
evaluation_metrics:
  - mean_absolute_error
  - root_mean_squared_error
  - correlation_coefficient
  - outlier_analysis
  
replication_requirements:
  minimum_participants: 3
  identical_environments: "enforced_via_containers"
  statistical_agreement: "MAE_within_0.1_across_replicators"
```

### **Sprint #2 (November 2025): DFT Conformer Energies**
```yaml
target_claim: "B3LYP/6-31G* reproduces experimental conformer ratios within 1 kcal/mol"
controversy: "basis_set_and_dispersion_correction_dependencies"
test_scope: "20_molecules_with_known_conformer_populations"
expected_participants: 12
computational_requirements: "gaussian_psi4_or_orca_access"
```

### **Sprint #3 (December 2025): Catalyst Screening**
```yaml
target_claim: "Machine learning models predict catalyst performance with R² > 0.8"
controversy: "dataset_leakage_and_generalization_concerns"
test_scope: "cross_validation_on_open_catalyst_datasets"
expected_participants: 20
computational_requirements: "python_ML_stack_gpu_recommended"
```

## 🏆 Truth Collapse Framework

### **Statistical Validation Criteria**
```python
class TruthCollapseValidator:
    def evaluate_claim(self, replications):
        """Determine if claim can be collapsed to truth."""
        
        criteria = {
            'statistical_power': len(replications) >= 3,
            'agreement_threshold': self.coefficient_of_variation(replications) < 0.15,
            'outlier_robustness': self.outlier_test(replications),
            'method_consistency': self.cross_method_agreement(replications),
            'environmental_reproducibility': self.environment_check(replications)
        }
        
        confidence_score = self.calculate_confidence(criteria)
        
        if confidence_score > 0.85:
            return self.collapse_to_truth(replications)
        elif confidence_score > 0.60:
            return self.mark_as_preliminary(replications)
        else:
            return self.flag_as_disputed(replications)
```

### **Truth Categories**
```yaml
validated_truth:
  confidence: ">85%"
  replications: "≥3 independent"
  agreement: "CV < 15%"
  badge: "🏆 Truth Validated"
  
preliminary_evidence:
  confidence: "60-85%"
  replications: "≥2 independent"
  agreement: "CV < 25%"
  badge: "⚠️ Preliminary Evidence"
  
disputed_claim:
  confidence: "<60%"
  replications: "conflicting results"
  agreement: "CV > 25%"
  badge: "❌ Disputed - Needs Investigation"
  
insufficient_data:
  confidence: "unknown"
  replications: "<2 independent"
  agreement: "unknown"
  badge: "❓ Insufficient Data"
```

## 📊 Community Recognition System

### **Contributor Levels**
```yaml
replication_badges:
  bronze_replicator: "1-5_successful_replications"
  silver_replicator: "6-15_successful_replications"  
  gold_replicator: "16+_successful_replications"
  
specialized_roles:
  method_expert: "deep_expertise_in_specific_computational_method"
  protocol_developer: "creates_standardized_testing_protocols"
  statistical_validator: "expert_in_replication_statistics"
  truth_adjudicator: "experienced_in_claim_evaluation"
  
leadership_recognition:
  sprint_organizer: "leads_monthly_sprint_events"
  community_champion: "exceptional_contribution_to_reproducibility"
  truth_pioneer: "first_to_collapse_disputed_claims"
```

### **Hall of Replications**
```markdown
## 🏆 Champion Replicators (Updated Monthly)

### October 2025 Sprint Champions
1. **Dr. Sarah Chen** - 5 successful replications, statistical validator
2. **Alex Rodriguez** - 4 replications, method expert (RDKit)
3. **Prof. Michael Kim** - 3 replications, protocol developer

### All-Time Truth Collapse Leaders
1. **Dr. Sarah Chen** - 12 collapsed claims, 95% replication success rate
2. **Dr. James Wilson** - 8 collapsed claims, statistical validator
3. **Prof. Lisa Garcia** - 7 collapsed claims, method expert
```

## 🔧 Technical Infrastructure

### **Locked Environments**
```yaml
container_registry: "ghcr.io/fotchemistry/reproducibility"

available_images:
  rdkit_environment:
    base: "ubuntu:20.04"
    python: "3.9.18"
    rdkit: "2023.09.1"
    packages: ["pandas", "numpy", "scipy", "matplotlib"]
    
  quantum_chemistry:
    base: "ubuntu:20.04"
    software: ["psi4:1.7", "xtb:6.6.1"]
    python: "3.9.18"
    packages: ["cclib", "ase", "qcengine"]
    
  machine_learning:
    base: "pytorch/pytorch:2.1.0-cuda11.8-cudnn8-runtime"
    packages: ["scikit-learn", "dgl", "torch-geometric"]
    models: ["pre_trained_catalyst_models"]
```

### **Results Database**
```yaml
database_schema:
  replications:
    - sprint_id
    - claim_id
    - replicator_id
    - method_used
    - environment_hash
    - results_data
    - timestamp
    - validation_status
    
  claims:
    - original_paper_doi
    - claim_statement
    - target_metrics
    - truth_status
    - confidence_score
    - collapse_timestamp
    
  participants:
    - user_id
    - orcid
    - institution
    - expertise_areas
    - reputation_score
```

## 📈 Sprint Analytics Dashboard

### **Real-Time Metrics** (Coming Soon)
Track progress at: [reproducibility.fotchemistry.org](https://reproducibility.fotchemistry.org)

```yaml
current_sprint_status:
  active_sprint: "Sprint #1 - pKa Predictions"
  days_remaining: 18
  registered_participants: 0
  completed_replications: 0
  preliminary_results: 0

all_time_statistics:
  total_sprints_completed: 0
  claims_tested: 0
  successful_replications: 0
  truth_collapses: 0
  disputed_resolutions: 0
  community_participants: 0
```

## 🚀 Participation Guide

### **For Individual Researchers**
1. **Register** for monthly sprint announcements
2. **Choose your expertise** area and preferred methods
3. **Join active sprint** when it matches your skills
4. **Follow protocols** exactly using provided containers
5. **Submit results** through standardized interface
6. **Earn recognition** through successful replications

### **For Research Groups**
1. **Organize team participation** in sprints
2. **Contribute computational resources** for complex claims
3. **Propose challenging claims** from your field
4. **Host sprint events** at your institution
5. **Mentor students** in reproducible research practices

### **For Institutions**
1. **Support faculty participation** in sprint activities
2. **Provide computational infrastructure** for complex replications
3. **Recognize reproducibility efforts** in promotion criteria
4. **Host community events** and workshops
5. **Contribute to truth collapse validation** processes

## 🎓 Educational Impact

### **Training Components**
- **Reproducibility Fundamentals**: Best practices for computational research
- **Statistical Validation**: Understanding replication statistics and confidence
- **Software Engineering**: Version control, containers, and environment management
- **Scientific Communication**: Reporting methods and results transparently

### **Student Opportunities**
- **Sprint Participation**: Hands-on experience with literature validation
- **Method Development**: Creating new replication protocols
- **Data Analysis**: Statistical evaluation of replication results
- **Community Leadership**: Organizing and leading sprint events

---

**Join the Reproducibility Sprint and help restore trust in the chemical literature through systematic validation and truth collapse.**

**Get Started**: [Register for Next Sprint](./register/) | [Nominate Claims](./nominate/) | [Access Environments](./environments/) | [View Hall of Fame](./hall_of_fame/)

*Monthly sprints organized by FoTChemistry - Open lab notebook & truth ledger for chemistry*
