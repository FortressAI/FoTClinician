#!/usr/bin/env python3
"""
FoTClinician Quick Accuracy Validation

Simplified validation test to demonstrate current system accuracy.
"""

import json
import sys
import time
import os

# Add project root to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from core.clinical.context_aware_quantum_engine import ContextAwareQuantumClinicalEngine

def run_quick_validation():
    """Run quick validation test"""
    
    print("🔬 FoTClinician Quick Accuracy Validation")
    print("=" * 50)
    
    # Initialize quantum engine
    engine = ContextAwareQuantumClinicalEngine(vqbit_dimension=256)
    
    # Test cases
    test_cases = [
        {
            "name": "Myocardial Infarction",
            "data": {
                "chief_complaint": "crushing chest pain radiating to left arm",
                "age": 55,
                "vital_signs": {"systolic_bp": 160, "heart_rate": 110},
                "symptoms": {"chest_pain": {"intensity": 1.0}}
            },
            "expected": "mi_acute"
        },
        {
            "name": "Diabetic Ketoacidosis", 
            "data": {
                "chief_complaint": "nausea, vomiting, confusion",
                "age": 72,
                "medical_history": ["diabetes_mellitus"],
                "laboratory": {"glucose": 485, "ph": 7.08}
            },
            "expected": "diabetic_ketoacidosis"
        },
        {
            "name": "Pediatric Sepsis",
            "data": {
                "chief_complaint": "fever, fussiness",
                "age": 3,  # 3 months
                "age_unit": "months",
                "vital_signs": {"temperature_c": 39.1, "heart_rate": 165},
                "laboratory": {"white_blood_count": 20.5}
            },
            "expected": "sepsis_infant"
        }
    ]
    
    results = []
    total_cases = len(test_cases)
    passed_cases = 0
    
    for i, case in enumerate(test_cases, 1):
        print(f"\n📋 Test Case {i}: {case['name']}")
        
        start_time = time.time()
        
        # Run quantum analysis
        quantum_case = engine.encode_clinical_case(case["data"])
        quantum_claim = engine.apply_virtue_supervision(quantum_case)
        
        end_time = time.time()
        response_time = end_time - start_time
        
        # Get top diagnosis from context-aware claim
        top_diagnosis = quantum_claim.primary_diagnosis
        top_confidence = quantum_claim.probability
        
        # Check if passed
        passed = top_diagnosis == case["expected"]
        if passed:
            passed_cases += 1
        
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"   {status} Expected: {case['expected']}")
        print(f"   Actual: {top_diagnosis} (confidence: {top_confidence:.1%})")
        print(f"   Age Context: {quantum_claim.age_context.age_band.value}")
        print(f"   ICD-10 Code: {quantum_claim.terminology_binding.icd10cm if quantum_claim.terminology_binding else 'N/A'}")
        print(f"   Response Time: {response_time:.2f}s")
        
        if quantum_claim.guidance_card:
            print(f"   ⚠️ Guidance Required: {quantum_claim.guidance_card.question}")
        
        results.append({
            "case_name": case["name"],
            "expected": case["expected"],
            "actual": top_diagnosis,
            "confidence": top_confidence,
            "passed": passed,
            "response_time": response_time,
            "age_context": quantum_claim.age_context.age_band.value,
            "icd10_code": quantum_claim.terminology_binding.icd10cm if quantum_claim.terminology_binding else "N/A",
            "guidance_required": quantum_claim.guidance_card is not None,
            "audit_trail": quantum_claim.audit_trail
        })
    
    # Summary
    accuracy = passed_cases / total_cases
    avg_response_time = sum(r["response_time"] for r in results) / len(results)
    avg_confidence = sum(r["confidence"] for r in results) / len(results)
    
    print(f"\n📊 VALIDATION SUMMARY")
    print(f"{'=' * 30}")
    print(f"🎯 Overall Accuracy: {accuracy:.1%} ({passed_cases}/{total_cases})")
    print(f"⚡ Average Response Time: {avg_response_time:.2f}s")
    print(f"🎓 Average Confidence: {avg_confidence:.1%}")
    
    # System status
    if accuracy >= 0.8:
        print(f"\n🎉 SYSTEM STATUS: EXCELLENT")
        print(f"   Clinical accuracy meets high standards")
    elif accuracy >= 0.6:
        print(f"\n✅ SYSTEM STATUS: GOOD")
        print(f"   Clinical accuracy is acceptable")
    else:
        print(f"\n⚠️ SYSTEM STATUS: NEEDS IMPROVEMENT")
        print(f"   Clinical accuracy below target")
    
    return {
        "accuracy": accuracy,
        "passed_cases": passed_cases,
        "total_cases": total_cases,
        "avg_response_time": avg_response_time,
        "avg_confidence": avg_confidence,
        "results": results
    }

if __name__ == "__main__":
    results = run_quick_validation()
    
    # Save results
    with open("quick_validation_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\n💾 Results saved to: quick_validation_results.json")
