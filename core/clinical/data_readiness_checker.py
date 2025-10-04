"""
FoT Clinical Data Readiness Checker

Validates clinical cases against data requirements for:
1. Medication Safety & Reconciliation  
2. Triage Assessment
3. Next Diagnostic Step (advice-level)

Follows Field of Truth principles:
- Every validation is a Claim with measurements and uncertainty
- NearMiss when data insufficient 
- Machine-checkable contracts with precise gaps flagged
"""

import json
from typing import Dict, List, Any, Tuple
from dataclasses import dataclass, asdict
from enum import Enum

class ValidationTrack(Enum):
    MEDICATION_SAFETY = "MedicationSafety"
    TRIAGE_ASSESSMENT = "TriageAssessment" 
    NEXT_DIAGNOSTIC_STEP = "NextDiagnosticStep"

class ValidationResult(Enum):
    READY = "READY"
    NOT_READY = "NOT_READY"
    NEAR_MISS = "NEAR_MISS"

@dataclass
class DataGap:
    field: str
    reason: str
    severity: str  # "required", "highly_recommended", "optional"
    example_value: str

@dataclass 
class ValidationClaim:
    """FoT Claim for data readiness validation"""
    measurement_type: str
    value: Any
    uncertainty: float
    provenance: Dict[str, str]
    collapse_policy: str
    timestamp: str
    toolchain_hash: str

@dataclass
class TrackValidationResult:
    track: ValidationTrack
    result: ValidationResult
    gaps: List[DataGap]
    success_threshold: float
    actual_score: float
    validation_claim: ValidationClaim
    ready_for_collapse: bool

class ClinicalDataContractValidator:
    """Validates clinical cases against FoT data contracts"""
    
    def __init__(self):
        self.validation_history = []
        self.toolchain_version = "FoT-Clinical-v1.0"
        
    def validate_case(self, case_data: Dict[str, Any]) -> List[TrackValidationResult]:
        """
        Validate clinical case across all tracks
        
        Returns List[TrackValidationResult] with detailed gaps and readiness scores
        """
        results = []
        
        # Validate each track
        results.append(self._validate_medication_safety(case_data))
        results.append(self._validate_triage_assessment(case_data))
        results.append(self._validate_next_diagnostic_step(case_data))
        
        return results
    
    def _validate_medication_safety(self, case: Dict[str, Any]) -> TrackValidationResult:
        """Validate medication safety track requirements"""
        gaps = []
        score = 0.0
        required_fields = 5  # medications, allergies, renal, hepatic, pregnancy
        
        # Check medications list
        if not self._has_valid_field(case, 'medications', lambda x: isinstance(x, list) and len(x) > 0):
            gaps.append(DataGap(
                field="medications",
                reason="Required list of current medications with name, dose, frequency",
                severity="required",
                example_value='[{"name": "metformin", "dose": "500mg", "frequency": "bid"}]'
            ))
        else:
            score += 1.0
            
        # Check allergies list (can be empty but must be present)
        if 'allergies' not in case:
            gaps.append(DataGap(
                field="allergies", 
                reason="Must explicitly state 'allergies': [] even if none",
                severity="required",
                example_value='{"allergies": []}'
            ))
        else:
            score += 1.0
            
        # Check renal function 
        if not self._has_valid_field(case, 'renal_function', self._valid_renal_function):
            gaps.append(DataGap(
                field="renal_function",
                reason="Required: eGFR (preferred) or creatinine within 12 months", 
                severity="required",
                example_value='{"egfr": 75.2, "unit": "mL/min/1.73m2"}'
            ))
        else:
            score += 1.0
            
        # Check hepatic panel
        if not self._has_valid_field(case, 'hepatic_panel', self._valid_hepatic_panel):
            gaps.append(DataGap(
                field="hepatic_panel",
                reason="Required: AST, ALT, bilirubin within 12 months",
                severity="required", 
                example_value='{"ast": 25, "alt": 30, "bilirubin": 1.2}'
            ))
        else:
            score += 1.0
            
        # Check pregnancy status (conditional)
        patient_age = case.get('age', 0)
        patient_sex = case.get('sex_at_birth', '').lower()
        
        if patient_sex == 'female' and 12 <= patient_age <= 55:
            if 'pregnancy_status' not in case:
                gaps.append(DataGap(
                    field="pregnancy_status",
                    reason="Required for female patients age 12-55",
                    severity="required",
                    example_value='{"pregnancy_status": "not_pregnant"}'
                ))
            else:
                score += 1.0
        else:
            score += 1.0  # Not applicable
            
        result = ValidationResult.READY if len(gaps) == 0 else ValidationResult.NOT_READY
        success_percentage = score / required_fields
        
        validation_claim = ValidationClaim(
            measurement_type="MedicationSafetyReadiness",
            value=success_percentage,
            uncertainty=0.05,  # 5% uncertainty in validation assessment
            provenance={
                "validator_version": self.toolchain_version,
                "contract_version": "2024.1",
                "validation_method": "field_presence_check"
            },
            collapse_policy="collapse_if_score >= 1.0",
            timestamp=self._timestamp(),
            toolchain_hash=self._calculate_hash(case)
        )
        
        return TrackValidationResult(
            track=ValidationTrack.MEDICATION_SAFETY,
            result=result,
            gaps=gaps,
            success_threshold=1.0,
            actual_score=success_percentage,
            validation_claim=validation_claim,
            ready_for_collapse=result == ValidationResult.READY
        )
    
    def _validate_triage_assessment(self, case: Dict[str, Any]) -> TrackValidationResult:
        """Validate triage assessment track requirements"""
        gaps = []
        score = 0.0
        required_fields = 6  # complaint + 5 vitals
        
        # Check chief complaint
        if not self._has_valid_field(case, 'chief_complaint', lambda x: isinstance(x, str) and len(x.strip()) > 0):
            gaps.append(DataGap(
                field="chief_complaint",
                reason="Required chief complaint description",
                severity="required",
                example_value='"Patient reports chest pain for 2 hours"'
            ))
        else:
            score += 1.0
            
        # Check vital signs
        vitals_required = ['systolic_bp', 'diastolic_bp', 'heart_rate', 'respiratory_rate', 'temperature_c', 'spo2']
        
        for vital in vitals_required:
            if not self._has_valid_field(case, vital, lambda x: isinstance(x, (int, float)) and x > 0):
                gaps.append(DataGap(
                    field=vital,
                    reason=f"Required vital sign measurement",
                    severity="required", 
                    example_value=f"{{'{vital}': 120}}" if vital != 'temperature_c' else "{'temperature_c': 37.1}"
                ))
            else:
                score += 1.0
                
        result = ValidationResult.READY if len(gaps) == 0 else ValidationResult.NOT_READY
        success_percentage = score / required_fields
        
        validation_claim = ValidationClaim(
            measurement_type="TriageAssessmentReadiness",
            value=success_percentage,
            uncertainty=0.05,
            provenance={
                "validator_version": self.toolchain_version,
                "contract_version": "2024.1", 
                "validation_method": "vital_signs_completeness"
            },
            collapse_policy="collapse_if_score >= 1.0 AND all_vitals_present",
            timestamp=self._timestamp(),
            toolchain_hash=self._calculate_hash(case)
        )
        
        return TrackValidationResult(
            track=ValidationTrack.TRIAGE_ASSESSMENT,
            result=result,
            gaps=gaps,
            success_threshold=1.0,
            actual_score=success_percentage,
            validation_claim=validation_claim,
            ready_for_collapse=result == ValidationResult.READY
        )
        
    def _validate_next_diagnostic_step(self, case: Dict[str, Any]) -> TrackValidationResult:
        """Validate next diagnostic step track requirements"""
        gaps = []
        score = 0.0
        required_fields = 2  # complaint + one of history/problems/labs
        
        # Check chief complaint
        if not self._has_valid_field(case, 'chief_complaint', lambda x: isinstance(x, str) and len(x.strip()) > 0):
            gaps.append(DataGap(
                field="chief_complaint",
                reason="Required chief complaint for diagnostic planning",
                severity="required",
                example_value='"Patient reports shortness of breath"'
            ))
        else:
            score += 1.0
            
        # Check for at least one: HPI, active problems, OR labs/imaging
        has_context = any([
            self._has_valid_field(case, 'hpi', lambda x: isinstance(x, str) and len(x.strip()) > 0),
            self._has_valid_field(case, 'active_problems', lambda x: isinstance(x, list) and len(x) > 0),
            self._has_valid_field(case, 'labs_imaging', lambda x: isinstance(x, str) and len(x.strip()) > 0)
        ])
        
        if not has_context:
            gaps.append(DataGap(
                field="clinical_context",
                reason="Need at least one: HPI, active problems, OR labs/imaging summary",
                severity="required",
                example_value='{"hpi": "Patient is 45M with chest pain x 2 hours, radiating to left arm"}'
            ))
        else:
            score += 1.0
            
        result = ValidationResult.READY if len(gaps) == 0 else ValidationResult.NOT_READY
        success_percentage = score / required_fields
        
        # Create NearMiss if ready but engines disagree
        if result == ValidationResult.READY:
            # Simulate engine disagreement check (in real system, would query actual engines)
            diagnostic_uncertainty = case.get('diagnostic_uncertainty', 0.3)
            if diagnostic_uncertainty > 0.25:  # High uncertainty threshold
                result = ValidationResult.NEAR_MISS
        
        validation_claim = ValidationClaim(
            measurement_type="NextDiagnosticReadiness", 
            value=success_percentage,
            uncertainty=diagnostic_uncertainty if 'diagnostic_uncertainty' in case else 0.1,
            provenance={
                "validator_version": self.toolchain_version,
                "contract_version": "2024.1",
                "validation_method": "clinical_context_sufficiency"
            },
            collapse_policy="collapse_if_score >= 1.0 AND uncertainty < 0.25", 
            timestamp=self._timestamp(),
            toolchain_hash=self._calculate_hash(case)
        )
        
        return TrackValidationResult(
            track=ValidationTrack.NEXT_DIAGNOSTIC_STEP,
            result=result,
            gaps=gaps,
            success_threshold=1.0,
            actual_score=success_percentage,
            validation_claim=validation_claim,
            ready_for_collapse=result == ValidationResult.READY
        )
    
    def _has_valid_field(self, case: Dict, field: str, validator_func) -> bool:
        """Check if field exists and passes validation function"""
        return field in case and validator_func(case[field])
    
    def _valid_renal_function(self, renal_data: Any) -> bool:
        """Validate renal function data"""
        if not isinstance(renal_data, dict):
            return False
        return 'egfr' in renal_data or ('creatinine' in renal_data and 'unit' in renal_data)
        
    def _valid_hepatic_panel(self, hepatic_data: Any) -> bool:
        """Validate hepatic panel data"""
        if not isinstance(hepatic_data, dict):
            return False
        return all(param in hepatic_data for param in ['ast', 'alt', 'bilirubin'])
        
    def generate_readiness_summary(self, validation_results: List[TrackValidationResult]) -> Dict[str, Any]:
        """Generate human-readable summary of validation results"""
        ready_tracks = [r.track.value for r in validation_results if r.result == ValidationResult.READY]
        not_ready_tracks = [r.track.value for r in validation_results if r.result == ValidationResult.NOT_READY]
        near_miss_tracks = [r.track.value for r in validation_results if r.result == ValidationResult.NEAR_MISS]
        
        total_gaps = sum(len(r.gaps) for r in validation_results)
        
        return {
            "overall_assessment": {
                "ready_tracks": ready_tracks,
                "not_ready_tracks": not_ready_tracks, 
                "near_miss_tracks": near_miss_tracks,
                "total_data_gaps": total_gaps
            },
            "track_details": [asdict(result) for result in validation_results],
            "recommendations": self._generate_recommendations(validation_results),
            "collapsed_claims": [r.validation_claim for r in validation_results if r.ready_for_collapse]
        }
    
    def _generate_recommendations(self, results: List[TrackValidationResult]) -> List[str]:
        """Generate specific recommendations for improving data readiness"""
        recommendations = []
        
        for result in results:
            if result.gaps:
                recommendations.append(f"TRACK {result.track.value}: Need {len(result.gaps)} more fields")
                for gap in result.gaps[:2]:  # Show top 2 gaps
                    recommendations.append(f"  - {gap.field}: {gap.reason}")
                    
        return recommendations
    
    def _timestamp(self) -> str:
        """Generate ISO timestamp"""
        from datetime import datetime
        return datetime.utcnow().isoformat() + "Z"
        
    def _calculate_hash(self, case_data: Dict) -> str:
        """Calculate hash for toolchain provenance"""
        import hashlib
        case_str = json.dumps(case_data, sort_keys=True)
        return hashlib.sha256(case_str.encode()).hexdigest()[:16]


def validate_clinical_case_from_file(case_file_path: str) -> None:
    """Command-line interface for validating clinical cases"""
    import argparse
    import sys
    
    try:
        with open(case_file_path, 'r') as f:
            case_data = json.load(f)
            
        validator = ClinicalDataContractValidator()
        results = validator.validate_case(case_data)
        summary = validator.generate_readiness_summary(results)
        
        print(f"\n🧬 FoT Clinical Data Readiness Validation")
        print(f"📁 Case: {case_file_path}")
        print(f"🕐 Validated: {summary['track_details'][0]['validation_claim']['timestamp']}")
        
        print(f"\n📊 Overall Assessment:")
        print(f"  ✅ Ready Tracks: {', '.join(summary['overall_assessment']['ready_tracks'])}")
        print(f"  ❌ Not Ready: {', '.join(summary['overall_assessment']['not_ready_tracks'])}")
        print(f"  ⚠️  Near Miss: {', '.join(summary['overall_assessment']['near_miss_tracks'])}")
        print(f"  📋 Total Data Gaps: {summary['overall_assessment']['total_data_gaps']}")
        
        if summary['recommendations']:
            print(f"\n💡 Recommendations:")
            for rec in summary['recommendations']:
                print(f"  {rec}")
                
        print(f"\n🔬 Collapsed Claims ({len(summary['collapsed_claims'])}):")
        for claim in summary['collapsed_claims']:
            print(f"  ✅ {claim['measurement_type']}: {claim['value']:.2f} ± {claim['uncertainty']:.3f}")
            
    except FileNotFoundError:
        print(f"❌ Case file not found: {case_file_path}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"❌ Invalid JSON in case file: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Validation error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python data_readiness_checker.py <case_file.json>")
        sys.exit(1)
        
    validate_clinical_case_from_file(sys.argv[1])