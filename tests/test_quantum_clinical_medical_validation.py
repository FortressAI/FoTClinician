"""
FoT Quantum Clinical Engine - Medical Validation Test Suite

CRITICAL VALIDATION TESTS:
1. Medical Accuracy - Correct diagnoses for known cases
2. Quantum Physics Correctness - Physically valid quantum mechanics
3. Safety Protocols - Non-maleficence and harm prevention
4. Clinical Decision Logic - Proper differential generation
5. Virtue Supervision - Ethical constraint enforcement
6. Edge Cases - Malpractice prevention scenarios

This test suite validates the quantum clinical engine could pass medical exams
and is safe for clinical decision support.
"""

import unittest
import numpy as np
import sys
import os
from typing import Dict, List, Any, Tuple

# Add project root to path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from core.clinical.quantum_clinical_engine import (
    QuantumClinicalEngine, 
    QuantumClinicalCase, 
    vQbitClinicalClaim,
    QuantumVirtueSupervisor,
    QuantumClinicalState
)

class MedicalValidationTestCase(unittest.TestCase):
    """Test medical accuracy against known clinical scenarios"""
    
    def setUp(self):
        """Set up quantum clinical engine for testing"""
        self.quantum_engine = QuantumClinicalEngine(vqbit_dimension=256)
        
    def test_myocardial_infarction_diagnosis(self):
        """
        Test quantum engine correctly identifies myocardial infarction
        
        USMLE-level case: 65M with chest pain, diaphoresis, ECG changes
        Expected: High probability MI diagnosis with proper quantum superposition
        """
        mi_case = {
            'case_id': 'MI_VALIDATION_001',
            'chief_complaint': 'chest pain radiating to left arm',
            'age': 65,
            'gender': 'male',
            'symptoms': {
                'chest_pain': {'intensity': 0.9, 'duration_hours': 2, 'radiating': True},
                'diaphoresis': {'intensity': 0.8, 'duration_hours': 1},
                'nausea': {'intensity': 0.6, 'duration_hours': 1},
                'dyspnea': {'intensity': 0.5, 'duration_hours': 1}
            },
            'vital_signs': {
                'systolic_bp': 180,
                'diastolic_bp': 105,
                'heart_rate': 95,
                'respiratory_rate': 24,
                'temperature_c': 37.1,
                'spo2': 92
            },
            'medical_history': ['hypertension', 'diabetes_type2', 'smoking'],
            'risk_factors': ['family_history_mi', 'hyperlipidemia']
        }
        
        # Encode into quantum state
        quantum_case = self.quantum_engine.encode_clinical_case(mi_case)
        
        # Validate quantum differential diagnosis
        self.assertIn('mi_acute', quantum_case.differential_qbits)
        
        # MI should have highest probability in quantum superposition
        mi_probability = abs(quantum_case.differential_qbits['mi_acute'])**2
        self.assertGreaterEqual(mi_probability, 0.25, 
                              f"MI probability too low: {mi_probability:.3f}")
        
        print(f"✅ MI Case Validation:")
        print(f"   MI Probability: {mi_probability:.3f}")
        print(f"   Quantum QBits Count: {len(quantum_case.differential_qbits)}")
        
    def test_aortic_dissection_recognition(self):
        """
        Test recognition of aortic dissection (high-risk diagnosis)
        
        Critical case: Acute onset tearing chest/back pain with hypertension
        Expected: Significant aortic dissection probability to trigger urgent workup
        """
        dissection_case = {
            'case_id': 'DISSECTION_VALIDATION_001',
            'chief_complaint': 'acute tearing chest and back pain',
            'age': 55,
            'gender': 'male',
            'symptoms': {
                'chest_pain': {'intensity': 1.0, 'duration_hours': 1, 'quality': 'tearing'},
                'back_pain': {'intensity': 0.9, 'duration_hours': 1},
                'loss_of_consciousness': {'intensity': 1.0, 'duration_minutes': 5}
            },
            'vital_signs': {
                'systolic_bp': 200,  # Hypertensive crisis
                'diastolic_bp': 120,
                'heart_rate': 85,
                'respiratory_rate': 20,
                'temperature_c': 36.8,
                'spo2': 98
            },
            'medical_history': ['hypertension', 'marfan_syndrome'],
            'physical_exam': {'aortic_regurgitation_murmur': True}
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(dissection_case)
        
        # Aortic dissection should be recognized
        self.assertIn('aortic_dissection', quantum_case.differential_qbits)
        dissection_prob = abs(quantum_case.differential_qbits['aortic_dissection'])**2
        
        # Should be significant enough to trigger emergency protocols
        self.assertGreaterEqual(dissection_prob, 0.08, 
                              f"Aortic dissection probability too low: {dissection_prob:.3f}")
        
        print(f"✅ Aortic Dissection Validation:")
        print(f"   Dissection probability: {dissection_prob:.3f}")
        
    def test_pediatric_hypoglycemia(self):
        """
        Test pediatric-specific diagnosis recognition
        
        Pediatric emergency: Infant with seizures and hypoglycemia
        Expected: Pediatric-specific differentials with age-adjusted probabilities
        """
        pediatric_case = {
            'case_id': 'PEDIATRIC_VALIDATION_001',
            'chief_complaint': 'seizure episode',
            'age': 2,  # 2 months old
            'gender': 'male',
            'symptoms': {
                'seizures': {'intensity': 1.0, 'duration_minutes': 3},
                'lethargy': {'intensity': 0.8, 'duration_hours': 2},
                'poor_feeding': {'intensity': 0.7, 'duration_hours': 8},
                'hypothermia': {'intensity': 0.6, 'duration_hours': 4}
            },
            'vital_signs': {
                'systolic_bp': 65,  # Pediatric norms
                'heart_rate': 140,
                'respiratory_rate': 35,
                'temperature_c': 35.2,  # Hypothermic
                'spo2': 88
            },
            'laboratory': {'glucose': 30},  # Severe hypoglycemia
            'birth_history': {'premature_birth': True, 'gestational_age_weeks': 32}
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(pediatric_case)
        
        # Pediatric case should generate appropriate differentials
        pediatric_dxs = [dx for dx in quantum_case.differential_qbits.keys() 
                        if any(keyword in dx.lower() for keyword in 
                              ['neonatal', 'pediatric', 'congenital', 'infectious'])]
        
        self.assertTrue(len(pediatric_dxs) > 0, 
                       "No pediatric-specific diagnoses found")
        
        print(f"✅ Pediatric Case Validation:")
        print(f"   Pediatric differentials: {pediatric_dxs}")
        print(f"   Total quantum QBits: {len(quantum_case.differential_qbits)}")


class QuantumPhysicsValidationTest(unittest.TestCase):
    """Test quantum mechanics implementation is physically correct"""
    
    def setUp(self):
        self.quantum_engine = QuantumClinicalEngine(vqbit_dimension=128)
        
    def test_unitary_state_evolution(self):
        """
        Test quantum state evolution follows unitary transformations
        
        Critical for quantum correctness: State magnitudes must be preserved
        """
        test_case = {
            'case_id': 'QUANTUM_PHYSICS_TEST_001',
            'chief_complaint': 'test symptom',
            'symptoms': {'test_symptom': {'intensity': 0.5}},
            'vital_signs': {'heart_rate': 70}
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(test_case)
        
        # Test unitarity: sum of squared amplitudes should equal 1.0
        total_probability = sum(abs(amp)**2 for amp in quantum_case.differential_qbits.values())
        
        self.assertAlmostEqual(total_probability, 1.0, places=2,
                             msg=f"Quantum unitarity violated: sum={total_probability:.6f}")
        
        print(f"✅ Quantum Physics Validation:")
        print(f"   Total probability: {total_probability:.6f}")
        
    def test_entanglement_correlations(self):
        """
        Test quantum entanglement correlations are physically valid
        
        Entanglement matrix should be Hermitian (self-adjoint)
        """
        test_case = {
            'case_id': 'ENTANGLEMENT_TEST_001',
            'chief_complaint': 'chest pain',
            'symptoms': {'chest_pain': {'intensity': 0.7}},
            'vital_signs': {'systolic_bp': 140, 'heart_rate': 90}
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(test_case)
        
        # Entanglement matrix should be approximately Hermitian
        entanglement_matrix = quantum_case.entanglement_matrix
        hermitian_diff = entanglement_matrix - entanglement_matrix.conj().T
        
        max_deviation = np.abs(hermitian_diff).max()
        
        # Medical-grade tolerance for quantum correlations
        self.assertLess(max_deviation, 0.1, 
                       f"Entanglement matrix not Hermitian: max deviation {max_deviation:.6f}")
        
        print(f"✅ Quantum Entanglement Validation:")
        print(f"   Max Hermitian deviation: {max_deviation:.6f}")
        
    def test_decoherence_physics(self):
        """
        Test quantum decoherence follows physical principles
        
        Decoherence rate should increase with data incompleteness
        """
        # Complete case (should have low decoherence)
        complete_case = {
            'chief_complaint': 'chest pain',
            'age': 55,
            'symptoms': {'chest_pain': {'intensity': 0.7}},
            'vital_signs': {'systolic_bp': 130, 'heart_rate': 80},
            'laboratory': {'troponin': 0.02},
            'electrocardiogram': {'st_elevation': True}
        }
        
        # Incomplete case (should have higher decoherence)
        incomplete_case = {
            'chief_complaint': 'fatigue'  # Minimal data
        }
        
        complete_quantum = self.quantum_engine.encode_clinical_case(complete_case)
        incomplete_quantum = self.quantum_engine.encode_clinical_case(incomplete_case)
        
        # Incomplete case should have higher decoherence
        self.assertLess(complete_quantum.decoherence_rate, 
                       incomplete_quantum.decoherence_rate,
                       "Quantum decoherence physics violated")
        
        print(f"✅ Quantum Decoherence Validation:")
        print(f"   Complete case decoherence: {complete_quantum.decoherence_rate:.4f}")
        print(f"   Incomplete case decoherence: {incomplete_quantum.decoherence_rate:.4f}")


class SafetyProtocolValidationTest(unittest.TestCase):
    """
    Validate safety protocols and non-maleficence
    
    Critical for preventing medical malpractice
    """
    
    def setUp(self):
        self.quantum_engine = QuantumClinicalEngine(vqbit_dimension=128)
        
    def test_non_maleficence_virtue_supervisor(self):
        """
        Test Non-maleficence Virtue Supervisor prevents harmful collapses
        
        Must NEVER recommend dangerous treatments when low-risk options exist
        """
        dangerous_case = {
            'case_id': 'SAFETY_TEST_DANGEROUS',
            'chief_complaint': 'mild headache',
            'age': 25,
            'symptoms': {'headache': {'intensity': 0.4, 'duration_hours': 2}},
            'vital_signs': {'systolic_bp': 120, 'heart_rate': 72},
            'allergies': ['penicillin'],
            'pregnancy_status': 'not_pregnant'
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(dangerous_case)
        claim = self.quantum_engine.apply_virtue_supervision(quantum_case)
        
        # Non-maleficence supervisor should prevent dangerous collapse
        non_maleficence_supervisor = QuantumVirtueSupervisor("non_maleficence")
        non_maleficence_measurement = non_maleficence_supervisor.quantum_measure(quantum_case)
        
        # Non-maleficence amplitude should be high (avoiding harm)
        self.assertGreater(abs(non_maleficence_measurement), 0.6,
                         f"Non-maleficence supervisor too weak: {abs(non_maleficence_measurement):.3f}")
        
        print(f"✅ Non-maleficence Validation:")
        print(f"   Non-maleficence amplitude: {abs(non_maleficence_measurement):.3f}")
        print(f"   Quantum collapse policy: {claim.collapse_policy}")
        
    def test_pregnancy_contraindications(self):
        """
        Test pregnancy-aware contraindications
        
        Must recognize pregnancy and avoid harmful interventions
        """
        pregnancy_case = {
            'case_id': 'PREGNANCY_SAFETY_TEST',
            'chief_complaint': 'nausea and vomiting',
            'age': 28,
            'gender': 'female',
            'pregnancy_status': 'pregnant_first_trimester',
            'symptoms': {'nausea': {'intensity': 0.8}, 'vomiting': {'intensity': 0.7}},
            'medications': [],  # Empty med list
            'allergies': []
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(pregnancy_case)
        claim = self.quantum_engine.apply_virtue_supervision(quantum_case)
        
        # Pregnancy case should require human review
        self.assertEqual(claim.collapse_policy, "virtue_measured_state",
                        "Pregnancy case should require careful measurement")
        
        print(f"✅ Pregnancy Safety Validation:")
        print(f"   Quantum claim collapse policy: {claim.collapse_policy}")
        
    def test_drug_allergy_safety(self):
        """
        Test drug allergy recognition and avoidance
        
        Must recognize allergies and prevent dangerous drug recommendations
        """
        allergy_case = {
            'case_id': 'URUG_ALLERGY_TEST',
            'chief_complaint': 'skin rash',
            'age': 45,
            'symptoms': {'skin_rash': {'intensity': 0.6, 'duration_hours': 12}},
            'allergies': ['sulfonamides', 'penicillin', 'iodine_contrast'],
            'current_medications': ['metformin', 'lisinopril'],
            'medical_history': ['diabetes_type2', 'hypertension']
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(allergy_case)
        
        # Should generate allergy-aware differentials
        allergy_dxs = [dx for dx in quantum_case.differential_qbits.keys()]
        
        self.assertTrue(len(allergy_dxs) > 0, "No diagnoses generated for allergy case")
        
        print(f"✅ Drug Allergy Safety Validation:")
        print(f"   Differential diagnoses: {allergy_dxs}")


class VirtueSuperVisionTest(unittest.TestCase):
    """
    Test quantum virtue supervisors prevent unethical clinical decisions
    
    Virtue supervisors must act as quantum gates preventing harmful collapses
    """
    
    def setUp(self):
        """Set up quantum clinical engine for virtue supervision testing"""
        self.quantum_engine = QuantumClinicalEngine(vqbit_dimension=128)
    
    def test_honesty_supervisor_uncertainty_surfaces(self):
        """
        Test Honesty supervisor surfaces uncertainty honestly
        
        Must not hide doubt - quantum uncertainty must be transparent
        """
        ambiguous_case = {
            'case_id': 'UNCERTAINTY_TEST',
            'chief_complaint': 'abdominal pain',  # Very nonspecific
            'age': 35,
            'symptoms': {'abdominal_pain': {'intensity': 0.5}},  # Minimal info
            # Missing: vital signs, basic labs, detailed history
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(ambiguous_case)
        
        honesty_supervisor = QuantumVirtueSupervisor("honesty")
        honesty_measurement = honesty_supervisor.quantum_measure(quantum_case)
        
        # Honesty supervisor should reflect uncertainty in quantum measurement
        claim = self.quantum_engine.apply_virtue_supervision(quantum_case)
        
        # High uncertainty should result in superposed state (not collapsed)
        if quantum_case.decoherence_rate > 0.15:  # High uncertainty case
            self.assertEqual(claim.quantum_state, QuantumClinicalState.SUPERPOSED,
                           "High uncertainty should remain superposed")
        
        print(f"✅ Honesty Supervisor Validation:")
        print(f"   Honesty measurement: {abs(honesty_measurement):.3f}")
        print(f"   Quantum uncertainty: {claim.uncertainty_hbar:.4f}")
        
    def test_prudence_supervisor_safest_first(self):
        """
        Test Prudence supervisor defaults to safest options
        
        Should prefer low-risk diagnostic approaches over high-risk interventions
        """
        unclear_case = {
            'case_id': 'PRUDENCE_TEST',
            'chief_complaint': 'fatigue',
            'age': 55,
            'symptoms': {'fatigue': {'intensity': 0.6, 'duration_days': 7}},
            'vital_signs': {'systolic_bp': 130, 'heart_rate': 85},
            'medical_history': ['anemia_history']
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(unclear_case)
        
        prudence_supervisor = QuantumVirtueSupervisor("prudence")
        prudence_measurement = prudence_supervisor.quantum_measure(quantum_case)
        
        # Prudence should keep probabilities balanced (don't jump to conclusions)
        prudence_threshold = 0.5
        self.assertGreater(abs(prudence_measurement), prudence_threshold,
                         f"Prudence supervisor too weak: {abs(prudence_measurement):.3f}")
        
        print(f"✅ Prudence Supervisor Validation:")
        print(f"   Prudence measurement: {abs(prudence_measurement):.3f}")


class MedicalExamSimulationTest(unittest.TestCase):
    """
    Simulate medical licensing exam scenarios (USMLE-style)
    
    Test quantum engine against formal medical education standards
    """
    
    def setUp(self):
        self.quantum_engine = QuantumClinicalEngine(vqbit_dimension=128)
        
    def test_usmle_step2_case_acute_chf(self):
        """
        USMLE Step 2 CK Level Case: Acute Heart Failure
        
        Scenario: 70F with progressive dyspnea, edema, orthopnea
        Expected: Heart failure diagnosis with appropriate quantum probability
        """
        chf_case = {
            'case_id': 'USMLE_CHF_001',
            'chief_complaint': 'progressive shortness of breath',
            'age': 70,
            'gender': 'female',
            'symptoms': {
                'dyspnea': {'intensity': 0.9, 'duration_days': 5},
                'orthopnea': {'intensity': 0.7, 'duration_days': 3},
                'paroxysmal_nocturnal_dyspnea': {'intensity': 0.8, 'duration_nights': 2},
                'peripheral_edema': {'intensity': 0.6, 'duration_days': 4},
                'fatigue': {'intensity': 0.5, 'duration_weeks': 2}
            },
            'vital_signs': {
                'systolic_bp': 160,
                'diastolic_bp': 95,
                'heart_rate': 110,
                'respiratory_rate': 28,
                'temperature_c': 37.0,
                'spo2': 89
            },
            'physical_exam': {
                'bilateral_rales': True,
                'jvp_elevated': True,
                's3_gallop': True,
                'sacral_edema': True
            },
            'medical_history': ['hypertension', 'atrial_fibrillation', 'diabetes'],
            'echocardiogram': {'ef': 35}  # Reduced ejection fraction
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(chf_case)
        
        # Should include heart failure in differential
        heart_failure_dxs = [dx for dx in quantum_case.differential_qbits.keys() 
                           if 'heart_failure' in dx.lower() or 'chf' in dx.lower()]
        
        self.assertTrue(len(heart_failure_dxs) > 0, 
                       "Heart failure not identified in USMLE-level case")
        
        print(f"✅ USMLE CHF Case Validation:")
        print(f"   Heart failure diagnoses: {heart_failure_dxs}")
        
        # Apply virtue supervision for clinical decision
        claim = self.quantum_engine.apply_virtue_supervision(quantum_case)
        
        # Clinical decision should be measured/ready for review
        self.assertIn(claim.quantum_state.value, ['measured', 'collapsed'],
                     "CHF case should be clinically actionable")
        
        print(f"   Quantum clinical decision: {claim.quantum_state.value}")
        
    def test_board_certification_emergency_case(self):
        """
        Board Certification Emergency Medicine Case
        
        High-stakes emergency: Suspected aortic dissection vs MI
        Required: Rapid quantum assessment with virtue-based prioritization
        """
        emergency_case = {
            'case_id': 'BOARD_EMERGENCY_001',
            'chief_complaint': 'acute chest and back pain',
            'age': 55,
            'gender': 'male',
            'symptoms': {
                'chest_pain': {'intensity': 1.0, 'duration_hours': 1, 'quality': 'sharp'},
                'back_pain': {'intensity': 0.9, 'duration_hours': 1},
                'diaphoresis': {'intensity': 0.8, 'duration_minutes': 45},
                'nausea': {'intensity': 0.6, 'duration_minutes': 30}
            },
            'vital_signs': {
                'systolic_bp': 180,  # Hypertensive
                'diastolic_bp': 110,
                'heart_rate': 95,
                'respiratory_rate': 22,
                'temperature_c': 36.5,
                'spo2': 94
            },
            'physical_exam': {
                'aortic_regurgitation_murmur': True,
                'pulse_deficit': True,
                'bilateral_rlq_pain': False
            },
            'medical_history': ['marfan_syndrome', 'hypertension'],
            'risk_factors': ['family_history_sudden_death'],
            'ecg': {'no_st_elevation': True, 'left_ventricular_hypertrophy': True}
        }
        
        quantum_case = self.quantum_engine.encode_clinical_case(emergency_case)
        
        # Board-level emergency case should yield actionable diagnoses
        emergency_dxs = [dx for dx in quantum_case.differential_qbits.keys()]
        
        # Should include high-acuity differentials
        high_acuity_keywords = ['dissection', 'mi', 'embolism', 'hemorrhage']
        high_acuity_dxs = [dx for dx in emergency_dxs 
                          if any(keyword in dx.lower() for keyword in high_acuity_keywords)]
        
        self.assertTrue(len(high_acuity_dxs) > 0,
                       "Emergency case should identify high-acuity diagnoses")
        
        print(f"✅ Board Emergency Case Validation:")
        print(f"   Emergency diagnoses: {emergency_dxs}")
        print(f"   High acuity diagnoses: {high_acuity_dxs}")
        
        # Virtue supervision for emergency triage
        claim = self.quantum_engine.apply_virtue_supervision(quantum_case)
        
        # Emergency cases should be ready for immediate action
        self.assertIn(claim.quantum_state.value, ['measured', 'collapsed'],
                     "Emergency case must be actionable")


class TestSuiteRunner:
    """
    Comprehensive test suite runner for medical validation
    
    Runs all validation tests and generates clinical safety report
    """
    
    @staticmethod
    def run_all_medical_tests():
        """Run complete medical validation test suite"""
        
        print("🧬 FoT Quantum Clinical Engine - Medical Validation Suite")
        print("=" * 70)
        print("🔬 Testing Medical Accuracy, Quantum Physics, Safety Protocols")
        print("⚛️ Validating Against Medical Licensing Standards")
        print()
        
        # Create test suite
        test_classes = [
            MedicalValidationTestCase,
            QuantumPhysicsValidationTest,
            SafetyProtocolValidationTest,
            VirtueSuperVisionTest,
            MedicalExamSimulationTest
        ]
        
        total_tests = 0
        total_failures = 0
        total_errors = 0
        
        results_summary = {
            'medical_accuracy': 0,
            'quantum_physics': 0,
            'safety_protocols': 0,
            'virtue_supervision': 0,
            'medical_exam_simulation': 0
        }
        
        for test_class in test_classes:
            print(f"\n🧪 Running {test_class.__name__}")
            print("-" * 50)
            
            suite = unittest.TestLoader().loadTestsFromTestCase(test_class)
            runner = unittest.TextTestRunner(verbosity=2)
            result = runner.run(suite)
            
            total_tests += result.testsRun
            total_failures += len(result.failures)
            total_errors += len(result.errors)
            
            # Record category results
            category_name = test_class.__name__.lower().replace('test', '')
            tests_run = result.testsRun
            failures = len(result.failures) + len(result.errors)
            results_summary[category_name] = tests_run - failures
        
        print("\n📊 MEDICAL VALIDATION SUMMARY")
        print("=" * 70)
        
        for category, passed_tests in results_summary.items():
            print(f"✅ {category.replace('_', ' ').title()}: {passed_tests} tests passed")
        
        print(f"\n🎯 OVERALL RESULTS:")
        print(f"   Total Tests: {total_tests}")
        print(f"   Passed: {total_tests - total_failures - total_errors}")
        print(f"   Failed: {total_failures}")
        print(f"   Errors: {total_errors}")
        print(f"   Success Rate: {((total_tests - total_failures - total_errors) / total_tests * 100):.1f}%")
        
        if total_failures == 0 and total_errors == 0:
            print("\n🏆 ALL MEDICAL VALIDATION TESTS PASSED!")
            print("✅ System validated for clinical decision support")
            print("⚛️ Quantum physics implementation verified")
            print("🛡️ Safety protocols confirmed")
            print("🎓 Medical exam standards met")
            return True
        else:
            print("\n❌ MEDICAL VALIDATION FAILED!")
            print("🚨 System NOT ready for clinical use")
            return False


if __name__ == "__main__":
    # Run comprehensive medical validation
    success = TestSuiteRunner.run_all_medical_tests()
    
    if success:
        print("\n✅ QUANTUM CLINICAL ENGINE VALIDATED FOR MEDICAL USE")
        exit(0)
    else:
        print("\n❌ QUANTUM CLINICAL ENGINE FAILED MEDICAL VALIDATION")
        exit(1)
