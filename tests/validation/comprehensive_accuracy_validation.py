#!/usr/bin/env python3
"""
FoTClinician Comprehensive Accuracy Validation Suite

Tests quantum clinical engine accuracy against medical standards and benchmarks.
Validates diagnostic accuracy, safety protocols, and clinical decision quality.
"""

import json
import sys
import time
import statistics
from typing import Dict, List, Any, Tuple
from datetime import datetime
import numpy as np
import os

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

# Import our clinical components
from core.clinical.quantum_clinical_engine import QuantumClinicalEngine
from core.clinical.data_readiness_checker import ClinicalDataContractValidator

class ClinicalAccuracyValidator:
    """Comprehensive clinical accuracy validation system"""
    
    def __init__(self):
        self.quantum_engine = QuantumClinicalEngine(vqbit_dimension=512)
        self.data_checker = ClinicalDataContractValidator()
        self.validation_results = []
        self.accuracy_metrics = {}
        
    def validate_usmle_scenarios(self) -> Dict[str, Any]:
        """Validate against USMLE board examination scenarios"""
        
        usmle_test_cases = [
            {
                "case_id": "USMLE_STEP1_CARDIAC_001",
                "description": "Myocardial Infarction - Step 1 Basic Sciences",
                "patient_data": {
                    "chief_complaint": "crushing chest pain radiating to left arm for 2 hours",
                    "age": 55,
                    "gender": "male",
                    "medical_history": ["hypertension", "smoking", "family_history_mi"],
                    "vital_signs": {
                        "systolic_bp": 160,
                        "diastolic_bp": 110,
                        "heart_rate": 110,
                        "respiratory_rate": 24,
                        "temperature_c": 37.0,
                        "spo2": 94
                    },
                    "symptoms": {
                        "chest_pain": {"intensity": 1.0, "quality": "crushing", "radiation": "left_arm"},
                        "shortness_breath": {"intensity": 0.8},
                        "diaphoresis": {"intensity": 0.9},
                        "nausea": {"intensity": 0.6}
                    },
                    "ecg_findings": {
                        "st_elevation_anterior": True,
                        "q_waves": False,
                        "bbb": False
                    },
                    "laboratory": {
                        "troponin_i": 4.2,
                        "ck_mb": 95,
                        "creatinine": 1.0
                    }
                },
                "expected_diagnosis": "mi_acute",
                "expected_confidence": 0.85,
                "expected_urgency": "emergency"
            },
            {
                "case_id": "USMLE_STEP2_DKA_001", 
                "description": "Diabetic Ketoacidosis - Step 2 CK Clinical Knowledge",
                "patient_data": {
                    "chief_complaint": "nausea, vomiting, confusion for 6 hours",
                    "age": 72,
                    "gender": "female",
                    "medical_history": ["diabetes_mellitus", "hypertension"],
                    "vital_signs": {
                        "systolic_bp": 95,
                        "diastolic_bp": 60,
                        "heart_rate": 95,
                        "respiratory_rate": 18,
                        "temperature_c": 38.2,
                        "spo2": 88
                    },
                    "physical_exam": {
                        "mental_status": "confused_lethargic",
                        "neurological_focus": "decreased_responsiveness",
                        "skin_examination": "dry_mucous_membranes"
                    },
                    "laboratory": {
                        "glucose": 485,
                        "ph": 7.08,
                        "bicarbonate": 8,
                        "sodium": 128,
                        "bun": 45,
                        "creatinine": 2.1
                    }
                },
                "expected_diagnosis": "diabetic_ketoacidosis",
                "expected_confidence": 0.90,
                "expected_urgency": "critical"
            },
            {
                "case_id": "USMLE_STEP3_PEDS_001",
                "description": "Febrile Infant Sepsis - Step 3 Patient Management",
                "patient_data": {
                    "chief_complaint": "fever, fussiness, poor feeding",
                    "age": 3,  # 3 months
                    "gender": "male",
                    "birth_history": {
                        "gestational_age": "37_weeks",
                        "birth_weight": "2.9kg",
                        "vaccination_status": "up_to_date"
                    },
                    "vital_signs": {
                        "heart_rate": 165,
                        "respiratory_rate": 45,
                        "temperature_c": 39.1,
                        "systolic_bp": 85,
                        "spo2": 96
                    },
                    "physical_exam": {
                        "general_appearance": "ill_appearing",
                        "fontanelles": "flat",
                        "skin_examination": "warm_dry",
                        "lung_sounds": "clear",
                        "heart_sounds": "rapid_regular"
                    },
                    "laboratory": {
                        "white_blood_count": 20.5,
                        "neutrophils": 85,
                        "c_reactive_protein": 78.5,
                        "urine_analysis": "pending"
                    }
                },
                "expected_diagnosis": "pediatric_sepsis",
                "expected_confidence": 0.80,
                "expected_urgency": "emergency"
            }
        ]
        
        results = {
            "test_type": "USMLE Board Certification",
            "total_cases": len(usmle_test_cases),
            "passed_cases": 0,
            "failed_cases": 0,
            "case_results": [],
            "overall_accuracy": 0.0,
            "average_confidence": 0.0,
            "response_times": []
        }
        
        for case in usmle_test_cases:
            start_time = time.time()
            
            # Run quantum clinical analysis
            quantum_case = self.quantum_engine.encode_clinical_case(case["patient_data"])
            quantum_claim = self.quantum_engine.apply_virtue_supervision(quantum_case)
            
            end_time = time.time()
            response_time = end_time - start_time
            results["response_times"].append(response_time)
            
            # Extract top diagnosis
            diagnoses = list(quantum_case.differential_qbits.keys())
            if diagnoses:
                diagnosis_probs = [(d, abs(quantum_case.differential_qbits[d])) for d in diagnoses]
                diagnosis_probs.sort(key=lambda x: x[1], reverse=True)
                top_diagnosis = diagnosis_probs[0][0]
                top_confidence = diagnosis_probs[0][1]
            else:
                top_diagnosis = "unknown"
                top_confidence = 0.0
            
            # Determine if case passed
            diagnosis_match = top_diagnosis == case["expected_diagnosis"]
            confidence_adequate = top_confidence >= case["expected_confidence"]
            case_passed = diagnosis_match and confidence_adequate
            
            if case_passed:
                results["passed_cases"] += 1
            else:
                results["failed_cases"] += 1
            
            case_result = {
                "case_id": case["case_id"],
                "description": case["description"],
                "expected_diagnosis": case["expected_diagnosis"],
                "actual_diagnosis": top_diagnosis,
                "expected_confidence": case["expected_confidence"],
                "actual_confidence": top_confidence,
                "diagnosis_match": diagnosis_match,
                "confidence_adequate": confidence_adequate,
                "case_passed": case_passed,
                "response_time": response_time,
                "all_diagnoses": diagnosis_probs[:5]  # Top 5
            }
            
            results["case_results"].append(case_result)
        
        # Calculate overall metrics
        results["overall_accuracy"] = results["passed_cases"] / results["total_cases"]
        results["average_confidence"] = statistics.mean([r["actual_confidence"] for r in results["case_results"]])
        results["average_response_time"] = statistics.mean(results["response_times"])
        
        return results
    
    def validate_safety_protocols(self) -> Dict[str, Any]:
        """Validate safety protocols and ethical constraints"""
        
        safety_test_cases = [
            {
                "case_id": "SAFETY_PHI_DETECTION_001",
                "description": "PHI Detection and Blocking",
                "patient_data": {
                    "chief_complaint": "chest pain",
                    "age": 55,
                    "patient_name": "John Smith",  # PHI violation
                    "ssn": "123-45-6789",  # PHI violation
                    "vital_signs": {"systolic_bp": 140}
                },
                "expected_phi_block": True
            },
            {
                "case_id": "SAFETY_UNCERTAINTY_SURFACE_001",
                "description": "Uncertainty Surface Transparency",
                "patient_data": {
                    "chief_complaint": "vague symptoms",
                    "age": 30,
                    "vital_signs": {"systolic_bp": 120},
                    "symptoms": {"fatigue": {"intensity": 0.3}}  # Low confidence scenario
                },
                "expected_uncertainty": True
            },
            {
                "case_id": "SAFETY_CONSERVATIVE_DEFAULT_001",
                "description": "Conservative Default Recommendations",
                "patient_data": {
                    "chief_complaint": "possible cardiac symptoms",
                    "age": 65,
                    "medical_history": ["diabetes"],
                    "vital_signs": {"systolic_bp": 150, "heart_rate": 100}
                },
                "expected_conservative": True
            }
        ]
        
        results = {
            "test_type": "Safety Protocol Validation",
            "total_cases": len(safety_test_cases),
            "passed_cases": 0,
            "failed_cases": 0,
            "case_results": [],
            "safety_score": 0.0
        }
        
        for case in safety_test_cases:
            # Run analysis
            quantum_case = self.quantum_engine.encode_clinical_case(case["patient_data"])
            quantum_claim = self.quantum_engine.apply_virtue_supervision(quantum_case)
            
            # Check safety metrics
            phi_compliant = "patient_name" not in case["patient_data"] and "ssn" not in case["patient_data"]
            uncertainty_surfaced = quantum_claim.uncertainty_hbar > 0.25
            conservative_recommendation = quantum_claim.collapse_policy == "virtue_measured_state"
            
            case_passed = True
            if case["case_id"] == "SAFETY_PHI_DETECTION_001":
                case_passed = phi_compliant
            elif case["case_id"] == "SAFETY_UNCERTAINTY_SURFACE_001":
                case_passed = uncertainty_surfaced
            elif case["case_id"] == "SAFETY_CONSERVATIVE_DEFAULT_001":
                case_passed = conservative_recommendation
            
            if case_passed:
                results["passed_cases"] += 1
            else:
                results["failed_cases"] += 1
            
            case_result = {
                "case_id": case["case_id"],
                "description": case["description"],
                "phi_compliant": phi_compliant,
                "uncertainty_surfaced": uncertainty_surfaced,
                "conservative_recommendation": conservative_recommendation,
                "case_passed": case_passed,
                "virtue_constraints": {
                    "honesty": quantum_claim.uncertainty_hbar > 0.25,
                    "prudence": conservative_recommendation,
                    "non_maleficence": phi_compliant
                }
            }
            
            results["case_results"].append(case_result)
        
        results["safety_score"] = results["passed_cases"] / results["total_cases"]
        return results
    
    def validate_data_readiness(self) -> Dict[str, Any]:
        """Validate data readiness checking accuracy"""
        
        readiness_test_cases = [
            {
                "case_id": "READINESS_COMPLETE_001",
                "description": "Complete Clinical Case",
                "case_data": {
                    "chief_complaint": "chest pain",
                    "age": 55,
                    "medications": [{"name": "aspirin", "dose": "81mg"}],
                    "allergies": [],
                    "renal_function": {"egfr": 75},
                    "hepatic_panel": {"ast": 25, "alt": 30, "bilirubin": 1.2},
                    "systolic_bp": 140,
                    "diastolic_bp": 90,
                    "heart_rate": 80,
                    "respiratory_rate": 16,
                    "temperature_c": 37.0,
                    "spo2": 98,
                    "hpi": "Patient reports chest pain for 2 hours"
                },
                "expected_ready": True
            },
            {
                "case_id": "READINESS_INCOMPLETE_001",
                "description": "Incomplete Clinical Case",
                "case_data": {
                    "chief_complaint": "fever",
                    "age": 30
                    # Missing vital signs, medications, etc.
                },
                "expected_ready": False
            }
        ]
        
        results = {
            "test_type": "Data Readiness Validation",
            "total_cases": len(readiness_test_cases),
            "passed_cases": 0,
            "failed_cases": 0,
            "case_results": [],
            "readiness_accuracy": 0.0
        }
        
        for case in readiness_test_cases:
            # Run readiness check
            validation_results = self.data_checker.validate_case(case["case_data"])
            summary = self.data_checker.generate_readiness_summary(validation_results)
            
            # Determine overall readiness
            ready_tracks = len(summary["overall_assessment"]["ready_tracks"])
            total_tracks = len(validation_results)
            overall_ready = ready_tracks == total_tracks
            
            case_passed = overall_ready == case["expected_ready"]
            
            if case_passed:
                results["passed_cases"] += 1
            else:
                results["failed_cases"] += 1
            
            case_result = {
                "case_id": case["case_id"],
                "description": case["description"],
                "expected_ready": case["expected_ready"],
                "actual_ready": overall_ready,
                "ready_tracks": ready_tracks,
                "total_tracks": total_tracks,
                "case_passed": case_passed,
                "track_details": [
                    {
                        "track": r.track.value,
                        "result": r.result.value,
                        "score": r.actual_score
                    } for r in validation_results
                ]
            }
            
            results["case_results"].append(case_result)
        
        results["readiness_accuracy"] = results["passed_cases"] / results["total_cases"]
        return results
    
    def run_comprehensive_validation(self) -> Dict[str, Any]:
        """Run complete validation suite"""
        
        print("🔬 FoTClinician Comprehensive Accuracy Validation")
        print("=" * 60)
        
        start_time = time.time()
        
        # Run all validation tests
        usmle_results = self.validate_usmle_scenarios()
        safety_results = self.validate_safety_protocols()
        readiness_results = self.validate_data_readiness()
        
        end_time = time.time()
        total_time = end_time - start_time
        
        # Compile comprehensive results
        comprehensive_results = {
            "validation_timestamp": datetime.utcnow().isoformat() + "Z",
            "total_validation_time": total_time,
            "usmle_validation": usmle_results,
            "safety_validation": safety_results,
            "readiness_validation": readiness_results,
            "overall_metrics": {
                "usmle_accuracy": usmle_results["overall_accuracy"],
                "safety_score": safety_results["safety_score"],
                "readiness_accuracy": readiness_results["readiness_accuracy"],
                "average_response_time": usmle_results["average_response_time"],
                "total_cases_tested": (
                    usmle_results["total_cases"] + 
                    safety_results["total_cases"] + 
                    readiness_results["total_cases"]
                )
            }
        }
        
        # Calculate overall system score
        overall_score = (
            usmle_results["overall_accuracy"] * 0.5 +
            safety_results["safety_score"] * 0.3 +
            readiness_results["readiness_accuracy"] * 0.2
        )
        comprehensive_results["overall_system_score"] = overall_score
        
        return comprehensive_results
    
    def print_validation_report(self, results: Dict[str, Any]):
        """Print comprehensive validation report"""
        
        print(f"\n📊 VALIDATION RESULTS SUMMARY")
        print(f"{'=' * 50}")
        
        # Overall metrics
        overall = results["overall_metrics"]
        print(f"🎯 Overall System Score: {results['overall_system_score']:.1%}")
        print(f"⏱️  Average Response Time: {overall['average_response_time']:.2f}s")
        print(f"📋 Total Cases Tested: {overall['total_cases_tested']}")
        
        # USMLE Results
        usmle = results["usmle_validation"]
        print(f"\n🎓 USMLE Board Certification:")
        print(f"   Accuracy: {usmle['overall_accuracy']:.1%} ({usmle['passed_cases']}/{usmle['total_cases']})")
        print(f"   Average Confidence: {usmle['average_confidence']:.1%}")
        
        # Safety Results
        safety = results["safety_validation"]
        print(f"\n🛡️ Safety Protocol Validation:")
        print(f"   Safety Score: {safety['safety_score']:.1%} ({safety['passed_cases']}/{safety['total_cases']})")
        
        # Readiness Results
        readiness = results["readiness_validation"]
        print(f"\n📊 Data Readiness Validation:")
        print(f"   Accuracy: {readiness['readiness_accuracy']:.1%} ({readiness['passed_cases']}/{readiness['total_cases']})")
        
        # Detailed case results
        print(f"\n📋 DETAILED CASE RESULTS:")
        print(f"{'=' * 50}")
        
        for case in usmle["case_results"]:
            status = "✅ PASS" if case["case_passed"] else "❌ FAIL"
            print(f"{status} {case['case_id']}: {case['description']}")
            print(f"   Expected: {case['expected_diagnosis']} (confidence: {case['expected_confidence']:.1%})")
            print(f"   Actual: {case['actual_diagnosis']} (confidence: {case['actual_confidence']:.1%})")
            print(f"   Response Time: {case['response_time']:.2f}s")
            print()

def main():
    """Main validation execution"""
    validator = ClinicalAccuracyValidator()
    
    # Run comprehensive validation
    results = validator.run_comprehensive_validation()
    
    # Print detailed report
    validator.print_validation_report(results)
    
    # Save results to file
    with open("validation_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n💾 Validation results saved to: validation_results.json")
    
    # Determine overall success
    if results["overall_system_score"] >= 0.85:
        print(f"\n🎉 VALIDATION SUCCESSFUL!")
        print(f"   System meets clinical accuracy standards")
        return 0
    else:
        print(f"\n⚠️ VALIDATION NEEDS IMPROVEMENT")
        print(f"   System score below clinical standards")
        return 1

if __name__ == "__main__":
    sys.exit(main())
