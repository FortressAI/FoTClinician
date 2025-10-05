#!/usr/bin/env python3
"""
FoT Context-Aware Quantum Clinical Engine

Implements context-aware clinical reasoning with age bands, terminology bindings,
and guidance cards for critical ambiguities only.
"""

import numpy as np
from typing import Dict, List, Any, Tuple, Optional
from dataclasses import dataclass
from enum import Enum
import hashlib
from datetime import datetime

class AgeBand(Enum):
    """Age band classifications for clinical context"""
    NEONATE = "neonate"      # ≤28 days
    INFANT = "infant"        # 29-365 days  
    CHILD = "child"          # 1-12 years
    ADOLESCENT = "adolescent" # 13-17 years
    ADULT = "adult"          # ≥18 years

class CareSetting(Enum):
    """Clinical care setting context"""
    EMERGENCY = "emergency"
    ICU = "icu"
    OUTPATIENT = "outpatient"
    INPATIENT = "inpatient"
    SURGERY = "surgery"

@dataclass
class AgeContext:
    """Demographic and clinical context"""
    age_days: int
    gestational_age_weeks: Optional[int] = None
    age_band: Optional[AgeBand] = None
    pregnancy_status: Optional[str] = None
    care_setting: Optional[CareSetting] = None
    device_context: Optional[str] = None
    
    def __post_init__(self):
        if self.age_band is None:
            self.age_band = self._determine_age_band()
    
    def _determine_age_band(self) -> AgeBand:
        """Determine age band from age in days"""
        if self.age_days <= 28:
            return AgeBand.NEONATE
        elif self.age_days <= 365:
            return AgeBand.INFANT
        elif self.age_days <= 4380:  # 12 years
            return AgeBand.CHILD
        elif self.age_days <= 6205:  # 17 years
            return AgeBand.ADOLESCENT
        else:
            return AgeBand.ADULT

@dataclass
class TerminologyBinding:
    """Clinical terminology binding"""
    snomed_id: str
    snomed_concept: str
    icd10cm: str
    cpt: Optional[str] = None
    hcpcs: Optional[str] = None
    binding_confidence: float = 1.0
    age_validity: List[AgeBand] = None
    care_setting_validity: List[CareSetting] = None
    
    def __post_init__(self):
        if self.age_validity is None:
            self.age_validity = [AgeBand.ADULT]  # Default to adult
        if self.care_setting_validity is None:
            self.care_setting_validity = list(CareSetting)

@dataclass
class ClinicalRule:
    """Clinical decision rule"""
    rule_id: str
    preconditions: Dict[str, Any]
    exclusions: List[str]
    priority: int
    source: str
    action: str

@dataclass
class GuidanceTrigger:
    """Guidance card trigger conditions"""
    material_impact: bool
    delta_score_threshold: float
    missing_field: Optional[str] = None
    safety_critical: bool = False

@dataclass
class GuidanceCard:
    """Guidance card for critical ambiguities"""
    primary_diagnosis: str
    rationale: str
    alternative_diagnoses: List[Tuple[str, float]]
    question: str
    material_impact: str
    safety_implications: str

class ContextAwareQuantumClinicalEngine:
    """Context-aware quantum clinical engine with proper terminology bindings"""
    
    def __init__(self, vqbit_dimension: int = 512):
        self.vqbit_dimension = vqbit_dimension
        self.age_context = None
        self.terminology_bindings = self._initialize_terminology_bindings()
        self.clinical_rules = self._initialize_clinical_rules()
        self.guidance_triggers = self._initialize_guidance_triggers()
    
    def _initialize_terminology_bindings(self) -> Dict[str, TerminologyBinding]:
        """Initialize SNOMED CT to ICD-10-CM terminology bindings"""
        return {
            # Sepsis bindings with age validity
            "sepsis_neonatal": TerminologyBinding(
                snomed_id="91302008",
                snomed_concept="Neonatal sepsis (disorder)",
                icd10cm="P36.9",
                age_validity=[AgeBand.NEONATE],
                binding_confidence=0.95
            ),
            "sepsis_infant": TerminologyBinding(
                snomed_id="91302008", 
                snomed_concept="Sepsis (disorder)",
                icd10cm="A41.9",
                age_validity=[AgeBand.INFANT, AgeBand.CHILD, AgeBand.ADOLESCENT, AgeBand.ADULT],
                binding_confidence=0.90
            ),
            "sepsis_staph_aureus": TerminologyBinding(
                snomed_id="91302008",
                snomed_concept="Sepsis (disorder)",
                icd10cm="A41.01",
                age_validity=[AgeBand.INFANT, AgeBand.CHILD, AgeBand.ADOLESCENT, AgeBand.ADULT],
                binding_confidence=0.95
            ),
            "sepsis_device_related": TerminologyBinding(
                snomed_id="91302008",
                snomed_concept="Sepsis (disorder)", 
                icd10cm="T80.211A",
                age_validity=[AgeBand.INFANT, AgeBand.CHILD, AgeBand.ADOLESCENT, AgeBand.ADULT],
                binding_confidence=0.90
            ),
            # Myocardial infarction bindings
            "mi_acute": TerminologyBinding(
                snomed_id="22298006",
                snomed_concept="Myocardial infarction (disorder)",
                icd10cm="I21.9",
                cpt="99291",
                age_validity=[AgeBand.ADOLESCENT, AgeBand.ADULT],
                binding_confidence=0.95
            ),
            # Diabetic ketoacidosis bindings
            "diabetic_ketoacidosis": TerminologyBinding(
                snomed_id="127012008",
                snomed_concept="Diabetic ketoacidosis (disorder)",
                icd10cm="E11.10",
                age_validity=[AgeBand.CHILD, AgeBand.ADOLESCENT, AgeBand.ADULT],
                binding_confidence=0.95
            )
        }
    
    def _initialize_clinical_rules(self) -> List[ClinicalRule]:
        """Initialize clinical decision rules"""
        return [
            ClinicalRule(
                rule_id="block_neonatal_codes_if_age_over_28d",
                preconditions={"age_days": ">28"},
                exclusions=["P36.x"],
                priority=1,
                source="ICD-10-CM Guidelines",
                action="exclude_neonatal_codes"
            ),
            ClinicalRule(
                rule_id="organism_specific_sepsis",
                preconditions={"lab_culture_organism": "present"},
                exclusions=[],
                priority=2,
                source="Clinical Guidelines",
                action="prefer_organism_specific_codes"
            ),
            ClinicalRule(
                rule_id="device_related_sepsis",
                preconditions={"has_vascular_device": True, "timing_hours": "<=48"},
                exclusions=[],
                priority=3,
                source="Infection Control Guidelines",
                action="prefer_device_related_codes"
            ),
            ClinicalRule(
                rule_id="unspecified_sepsis_default",
                preconditions={"no_organism": True, "not_device_related": True},
                exclusions=[],
                priority=4,
                source="Default Coding",
                action="prefer_unspecified_sepsis"
            )
        ]
    
    def _initialize_guidance_triggers(self) -> List[GuidanceTrigger]:
        """Initialize guidance card triggers"""
        return [
            GuidanceTrigger(
                material_impact=True,
                delta_score_threshold=0.03,
                safety_critical=True,
                missing_field="device_related_status"
            ),
            GuidanceTrigger(
                material_impact=True,
                delta_score_threshold=0.05,
                safety_critical=False,
                missing_field="organism_identification"
            )
        ]
    
    def encode_clinical_case(self, clinical_data: Dict[str, Any]) -> 'ContextAwareQuantumCase':
        """Encode clinical case with context awareness"""
        
        # Extract and normalize context
        age_days = self._calculate_age_days(clinical_data)
        self.age_context = AgeContext(
            age_days=age_days,
            gestational_age_weeks=clinical_data.get('gestational_age_weeks'),
            pregnancy_status=clinical_data.get('pregnancy_status'),
            care_setting=self._determine_care_setting(clinical_data),
            device_context=clinical_data.get('device_context')
        )
        
        # Generate quantum state vector
        quantum_state_vector = np.random.normal(0, 1, self.vqbit_dimension) + 1j * np.random.normal(0, 1, self.vqbit_dimension)
        quantum_state_vector = quantum_state_vector / np.linalg.norm(quantum_state_vector)
        
        # Generate context-aware differential diagnoses
        differential_qbits = self._generate_context_aware_differentials(clinical_data)
        
        # Create entanglement matrix
        entanglement_matrix = np.random.normal(0, 0.1, (len(differential_qbits), len(differential_qbits)))
        
        return ContextAwareQuantumCase(
            case_id=hashlib.md5(str(clinical_data).encode()).hexdigest()[:16],
            quantum_state_vector=quantum_state_vector,
            differential_qbits=differential_qbits,
            entanglement_matrix=entanglement_matrix,
            age_context=self.age_context,
            clinical_data=clinical_data
        )
    
    def _calculate_age_days(self, clinical_data: Dict[str, Any]) -> int:
        """Calculate age in days from clinical data"""
        age = clinical_data.get('age', 50)
        
        # If age is in months (for infants), convert to days
        if age < 2 and 'months' in str(clinical_data.get('age_unit', '')).lower():
            return age * 30
        elif age < 18 and age > 1:  # Likely years
            return age * 365
        else:
            return age * 365  # Default to years
    
    def _determine_care_setting(self, clinical_data: Dict[str, Any]) -> CareSetting:
        """Determine care setting from clinical data"""
        chief_complaint = clinical_data.get('chief_complaint', '').lower()
        
        if any(term in chief_complaint for term in ['emergency', 'chest pain', 'severe']):
            return CareSetting.EMERGENCY
        elif any(term in chief_complaint for term in ['icu', 'critical', 'ventilator']):
            return CareSetting.ICU
        elif any(term in chief_complaint for term in ['surgery', 'surgical', 'operation']):
            return CareSetting.SURGERY
        else:
            return CareSetting.OUTPATIENT
    
    def _generate_context_aware_differentials(self, clinical_data: Dict[str, Any]) -> Dict[str, complex]:
        """Generate context-aware differential diagnoses with terminology bindings"""
        
        chief_complaint = clinical_data.get('chief_complaint', '').lower()
        age_band = self.age_context.age_band
        laboratory = clinical_data.get('laboratory', {})
        
        differential_qbits = {}
        
        # Apply clinical rules and age validity constraints
        if 'fever' in chief_complaint or 'sepsis' in chief_complaint:
            # Check age validity for sepsis codes
            if age_band == AgeBand.NEONATE:
                # Neonatal sepsis - high confidence for classic presentation
                differential_qbits['sepsis_neonatal'] = 0.85 + 0.1j
            else:
                # Non-neonatal sepsis - check for organism specificity
                if laboratory.get('white_blood_count', 0) > 15:
                    if laboratory.get('culture_organism'):
                        differential_qbits['sepsis_staph_aureus'] = 0.90 + 0.05j
                    else:
                        differential_qbits['sepsis_infant'] = 0.80 + 0.15j
                else:
                    differential_qbits['sepsis_infant'] = 0.75 + 0.20j
        
        elif 'chest pain' in chief_complaint:
            # Myocardial infarction - age validity check
            if age_band in [AgeBand.ADOLESCENT, AgeBand.ADULT]:
                # High confidence for classic MI presentation
                differential_qbits['mi_acute'] = 0.85 + 0.10j
            else:
                # Pediatric chest pain - different differentials
                differential_qbits['pediatric_chest_pain'] = 0.70 + 0.25j
        
        elif any(term in chief_complaint for term in ['nausea', 'vomiting', 'confusion']):
            # Diabetic ketoacidosis - age validity check
            if age_band in [AgeBand.CHILD, AgeBand.ADOLESCENT, AgeBand.ADULT]:
                if laboratory.get('glucose', 0) > 300:
                    # Very high confidence for classic DKA with hyperglycemia
                    differential_qbits['diabetic_ketoacidosis'] = 0.95 + 0.03j
                else:
                    differential_qbits['diabetic_ketoacidosis'] = 0.80 + 0.15j
        
        return differential_qbits
    
    def apply_virtue_supervision(self, quantum_case: 'ContextAwareQuantumCase') -> 'ContextAwareQuantumClaim':
        """Apply virtue supervision with context awareness"""
        
        # Get top diagnosis with terminology binding
        diagnoses = list(quantum_case.differential_qbits.keys())
        if diagnoses:
            diagnosis_probs = [(d, abs(quantum_case.differential_qbits[d])) for d in diagnoses]
            diagnosis_probs.sort(key=lambda x: x[1], reverse=True)
            top_diagnosis = diagnosis_probs[0][0]
            base_confidence = diagnosis_probs[0][1]
            
            # Apply confidence boosting based on clinical factors
            top_confidence = self._boost_confidence(base_confidence, top_diagnosis, quantum_case)
        else:
            top_diagnosis = "unknown"
            top_confidence = 0.0
        
        # Get terminology binding
        terminology_binding = self.terminology_bindings.get(top_diagnosis)
        
        # Check for guidance triggers
        guidance_card = self._check_guidance_triggers(quantum_case, diagnosis_probs)
        
        # Calculate uncertainty
        uncertainty_hbar = 1.0 - top_confidence
        
        # Determine collapse policy
        if guidance_card:
            collapse_policy = "guidance_required"
        elif uncertainty_hbar > 0.25:
            collapse_policy = "virtue_measured_state"
        else:
            collapse_policy = "quantum_collapsed"
        
        return ContextAwareQuantumClaim(
            measurement_type="context_aware_diagnosis",
            quantum_state="collapsed" if collapse_policy == "quantum_collapsed" else "superposed",
            amplitude=complex(top_confidence, uncertainty_hbar),
            probability=top_confidence,
            phase=0.0,
            entanglement_list=[],
            collapse_policy=collapse_policy,
            uncertainty_hbar=uncertainty_hbar,
            toolchain_hash=hashlib.md5(f"{top_diagnosis}{top_confidence}".encode()).hexdigest()[:16],
            timestamp=datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
            primary_diagnosis=top_diagnosis,
            terminology_binding=terminology_binding,
            age_context=quantum_case.age_context,
            guidance_card=guidance_card,
            audit_trail={
                "age_band": quantum_case.age_context.age_band.value,
                "care_setting": quantum_case.age_context.care_setting.value,
                "rules_applied": self._get_applied_rules(quantum_case),
                "excluded_diagnoses": self._get_excluded_diagnoses(quantum_case),
                "confidence_score": top_confidence,
                "uncertainty_level": uncertainty_hbar
            }
        )
    
    def _check_guidance_triggers(self, quantum_case: 'ContextAwareQuantumCase', 
                               diagnosis_probs: List[Tuple[str, float]]) -> Optional[GuidanceCard]:
        """Check if guidance card should be triggered"""
        
        if len(diagnosis_probs) < 2:
            return None
        
        top_score = diagnosis_probs[0][1]
        second_score = diagnosis_probs[1][1]
        delta = top_score - second_score
        
        # Check guidance triggers
        for trigger in self.guidance_triggers:
            if delta < trigger.delta_score_threshold and trigger.material_impact:
                
                primary_diagnosis = diagnosis_probs[0][0]
                alternative_diagnoses = diagnosis_probs[1:3]  # Top 2 alternatives
                
                # Generate guidance card
                if trigger.missing_field == "device_related_status":
                    return GuidanceCard(
                        primary_diagnosis=primary_diagnosis,
                        rationale=f"Age-appropriate diagnosis for {quantum_case.age_context.age_band.value}",
                        alternative_diagnoses=alternative_diagnoses,
                        question="Device-related bloodstream infection suspected?",
                        material_impact="Affects billing code (A41.9 vs T80.211A)",
                        safety_implications="Device-related infections require different treatment protocols"
                    )
                elif trigger.missing_field == "organism_identification":
                    return GuidanceCard(
                        primary_diagnosis=primary_diagnosis,
                        rationale=f"Sepsis diagnosis with unclear organism",
                        alternative_diagnoses=alternative_diagnoses,
                        question="Blood culture results available?",
                        material_impact="Affects organism-specific coding",
                        safety_implications="Organism-specific treatment may be required"
                    )
        
        return None
    
    def _get_applied_rules(self, quantum_case: 'ContextAwareQuantumCase') -> List[str]:
        """Get list of applied clinical rules"""
        applied_rules = []
        
        for rule in self.clinical_rules:
            if self._rule_applies(rule, quantum_case):
                applied_rules.append(rule.rule_id)
        
        return applied_rules
    
    def _rule_applies(self, rule: ClinicalRule, quantum_case: 'ContextAwareQuantumCase') -> bool:
        """Check if a clinical rule applies to the case"""
        # Simplified rule application logic
        if rule.rule_id == "block_neonatal_codes_if_age_over_28d":
            return quantum_case.age_context.age_days > 28
        
        return False
    
    def _get_excluded_diagnoses(self, quantum_case: 'ContextAwareQuantumCase') -> List[str]:
        """Get list of excluded diagnoses with reasons"""
        excluded = []
        
        # Age-based exclusions
        if quantum_case.age_context.age_days > 28:
            excluded.append("neonatal_sepsis (age > 28 days)")
        
        if quantum_case.age_context.age_band in [AgeBand.NEONATE, AgeBand.INFANT, AgeBand.CHILD]:
            excluded.append("mi_acute (age < 13 years)")
        
        return excluded
    
    def _boost_confidence(self, base_confidence: float, diagnosis: str, quantum_case: 'ContextAwareQuantumCase') -> float:
        """Boost confidence based on clinical factors"""
        
        boosted_confidence = base_confidence
        clinical_data = quantum_case.clinical_data
        vital_signs = clinical_data.get('vital_signs', {})
        laboratory = clinical_data.get('laboratory', {})
        
        # Age appropriateness boost
        if quantum_case.age_context.age_band in [AgeBand.ADULT, AgeBand.ADOLESCENT]:
            boosted_confidence += 0.05  # Adult presentations are more predictable
        
        # Vital signs correlation boost
        if diagnosis == 'mi_acute':
            if vital_signs.get('systolic_bp', 0) > 140:  # Hypertension supports MI
                boosted_confidence += 0.08
            if vital_signs.get('heart_rate', 0) > 100:  # Tachycardia supports MI
                boosted_confidence += 0.05
        
        elif diagnosis == 'diabetic_ketoacidosis':
            if laboratory.get('glucose', 0) > 400:  # Severe hyperglycemia
                boosted_confidence += 0.10
            if laboratory.get('ph', 0) < 7.2:  # Acidosis
                boosted_confidence += 0.08
        
        elif diagnosis == 'sepsis_infant':
            if vital_signs.get('temperature_c', 0) > 38.5:  # High fever
                boosted_confidence += 0.08
            if laboratory.get('white_blood_count', 0) > 20:  # High WBC
                boosted_confidence += 0.10
            if vital_signs.get('heart_rate', 0) > 160:  # Tachycardia
                boosted_confidence += 0.05
        
        # Terminology binding confidence boost
        terminology_binding = self.terminology_bindings.get(diagnosis)
        if terminology_binding and terminology_binding.binding_confidence > 0.9:
            boosted_confidence += 0.05
        
        # Cap confidence at 0.95 to maintain some uncertainty
        return min(boosted_confidence, 0.95)

@dataclass
class ContextAwareQuantumCase:
    """Context-aware quantum clinical case"""
    case_id: str
    quantum_state_vector: np.ndarray
    differential_qbits: Dict[str, complex]
    entanglement_matrix: np.ndarray
    age_context: AgeContext
    clinical_data: Dict[str, Any]

@dataclass
class ContextAwareQuantumClaim:
    """Context-aware quantum clinical claim with terminology binding"""
    measurement_type: str
    quantum_state: str
    amplitude: complex
    probability: float
    phase: float
    entanglement_list: List[str]
    collapse_policy: str
    uncertainty_hbar: float
    toolchain_hash: str
    timestamp: str
    primary_diagnosis: str
    terminology_binding: Optional[TerminologyBinding]
    age_context: AgeContext
    guidance_card: Optional[GuidanceCard]
    audit_trail: Dict[str, Any]
