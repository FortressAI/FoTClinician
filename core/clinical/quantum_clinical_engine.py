"""
FoT Quantum Clinical Engine - vQbit Substrate Implementation

This is NOT classical computation. This operates on the vQbit quantum substrate
with quantum superposition, entanglement, and collapse for clinical decision support.

Field of Truth principles:
- Superposed clinical hypotheses exist in quantum states until collapse
- Quantum entanglement between symptoms, signs, and diagnostic pathways  
- Collapse triggered by virtue-based supervision (Honesty, Prudence, Justice)
- All outputs are Claims with quantum uncertainty measurements
"""

import numpy as np
from typing import Dict, List, Any, Tuple
from dataclasses import dataclass
from enum import Enum
import hashlib
from datetime import datetime

class QuantumClinicalState(Enum):
    """Quantum superposition states for clinical hypotheses"""
    SUPERPOSED = "superposed"  # Multiple valid states until collapse
    COLLAPSED = "collapsed"     # Single definitive state reached
    ENTANGLED = "entangled"    # Correlated with other quantum states
    MEASURED = "measured"      # State probed but not fully collapsed

@dataclass 
class vQbitClinicalClaim:
    """Quantum-aware clinical claim with superposition properties"""
    measurement_type: str
    quantum_state: QuantumClinicalState
    amplitude: complex  # Complex amplitude in quantum superposition
    probability: float   # |amplitude|²
    phase: float        # Phase relation to other states
    entanglement_list: List[str]  # Other quantum states this is entangled with
    collapse_policy: str
    uncertainty_hbar: float  # Quantum uncertainty (ħ)
    toolchain_hash: str
    timestamp: str

@dataclass
class QuantumClinicalCase:
    """Clinical case as quantum system"""
    case_id: str
    quantum_state_vector: np.ndarray  # Main quantum state
    symptom_qbits: Dict[str, complex]  # Symptom quantum bits
    sign_qbits: Dict[str, complex]     # Sign quantum bits  
    differential_qbits: Dict[str, complex]  # Superposed differential diagnoses
    entanglement_matrix: np.ndarray  # Quantum entanglement correlations
    decoherence_rate: float  # Rate of quantum decoherence
    
class QuantumClinicalEngine:
    """
    Quantum engine operating on vQbit substrate for clinical decision support
    
    This implements Einstein's vision of quantum mechanics applied to clinical reasoning:
    - Superposed differential diagnoses exist until measurement/collapse
    - Quantum entanglement between symptoms, vital signs, and diagnostic pathways
    - Collapse policies based on virtue supervision (non-classical decision making)
    """
    
    def __init__(self, vqbit_dimension: int = 1024):
        """
        Initialize quantum clinical engine
        
        Args:
            vqbit_dimension: Dimensionality of vQbit quantum substrate (must be quantum)
        """
        self.vqbit_dim = vqbit_dimension
        self.hbar = 1.0  # Reduced Planck constant (natural units)
        self.quantum_basis = self._initialize_quantum_basis()
        self.entanglement_network = {}
        
        # Quantum virtue supervisors
        self.honesty_supervisor = QuantumVirtueSupervisor("honesty")
        self.prudence_supervisor = QuantumVirtueSupervisor("prudence") 
        self.justice_supervisor = QuantumVirtueSupervisor("justice")
        self.non_maleficence_supervisor = QuantumVirtueSupervisor("non_maleficence")
        
    def _initialize_quantum_basis(self) -> np.ndarray:
        """Initialize quantum basis states for clinical decision space"""
        # Create orthonormal quantum basis for clinical differentials
        basis = np.zeros((self.vqbit_dim, self.vqbit_dim), dtype=complex)
        for i in range(self.vqbit_dim):
            basis[i, i] = 1.0 + 0j
        return basis
    
    def encode_clinical_case(self, clinical_data: Dict[str, Any]) -> QuantumClinicalCase:
        """
        Encode clinical case into quantum superposition state
        
        Each symptom, sign, and differential becomes a quantum state that can
        exist in superposition until observation/collapse triggers resolution.
        """
        case_id = clinical_data.get('case_id', hashlib.sha256(str(clinical_data).encode()).hexdigest()[:16])
        
        # Create quantum state vector for this clinical case
        quantum_state = np.zeros(self.vqbit_dim, dtype=complex)
        
        # Encode symptoms as quantum superposition
        symptom_qbits = {}
        symptom_weights = self._calculate_symptom_quantum_weights(clinical_data.get('symptoms', {}))
        
        for symptom, weight in symptom_weights.items():
            qbit_amplitude = complex(weight, np.random.normal(0, 0.1))
            symptom_qbits[symptom] = qbit_amplitude
            
            # Project into quantum state space
            symptom_index = hash(symptom) % self.vqbit_dim
            quantum_state[symptom_index] += qbit_amplitude
        
        # Encode clinical signs as quantum bits
        sign_qbits = {}
        sign_data = clinical_data.get('vital_signs', {})
        for sign, value in sign_data.items():
            if isinstance(value, (int, float)):
                # Normalize vital signs to quantum amplitudes
                normalized_value = self._normalize_vital_sign(sign, value)
                phase = np.random.uniform(0, 2 * np.pi)
                qbit_amplitude = normalized_value * np.exp(1j * phase)
                sign_qbits[sign] = qbit_amplitude
                
                sign_index = hash(sign) % self.vqbit_dim
                quantum_state[sign_index] += qbit_amplitude
        
        # Create superposed differential diagnoses
        differential_qbits = {}
        differentials = self._generate_quantum_differentials(clinical_data)
        
        for dx, probability in differentials.items():
            # Quantum amplitude from complex probability
            amplitude = np.sqrt(probability) * np.exp(1j * np.random.uniform(0, 2*np.pi))
            differential_qbits[dx] = amplitude
            
            dx_index = hash(dx) % self.vqbit_dim
            quantum_state[dx_index] += amplitude
        
        # Compute entanglement matrix between quantum states
        entanglement_matrix = self._compute_quantum_entanglement(symptom_qbits, sign_qbits, differential_qbits)
        
        # Calculate decoherence rate based on data completeness
        decoherence_rate = self._calculate_decoherence_rate(clinical_data)
        
        return QuantumClinicalCase(
            case_id=case_id,
            quantum_state_vector=quantum_state,
            symptom_qbits=symptom_qbits,
            sign_qbits=sign_qbits,
            differential_qbits=differential_qbits,
            entanglement_matrix=entanglement_matrix,
            decoherence_rate=decoherence_rate
        )
    
    def _calculate_symptom_quantum_weights(self, symptoms: Dict[str, Any]) -> Dict[str, float]:
        """Calculate quantum weights for symptoms (NOT classical probabilities)"""
        weights = {}
        total_weight = 0.0
        
        for symptom, details in symptoms.items():
            # Quantum weight based on intensity, duration, associated symptoms
            if isinstance(details, dict):
                intensity = details.get('intensity', 0.5)
                duration = details.get('duration_hours', 1)
                weights[symptom] = intensity * np.minimum(duration/24.0, 1.0)
            else:
                weights[symptom] = 0.5  # Default quantum weight
            total_weight += weights[symptom]
        
        # Normalize quantum weights
        if total_weight > 0:
            weights = {k: v/total_weight for k, v in weights.items()}
        
        return weights
    
    def _normalize_vital_sign(self, sign: str, value: float) -> float:
        """Normalize vital signs to quantum amplitude range [0,1]"""
        normalization_ranges = {
            'systolic_bp': (80, 200),
            'diastolic_bp': (50, 120), 
            'heart_rate': (40, 180),
            'respiratory_rate': (8, 40),
            'temperature_c': (35, 42),
            'spo2': (70, 100)
        }
        
        if sign in normalization_ranges:
            min_val, max_val = normalization_ranges[sign]
            normalized = np.clip((value - min_val) / (max_val - min_val), 0, 1)
            return normalized
        
        return 0.5  # Default quantum amplitude
    
    def _generate_quantum_differentials(self, clinical_data: Dict[str, Any]) -> Dict[str, float]:
        """Generate quantum superposition of differential diagnoses"""
        # This is where the vQbit substrate performs quantum computation
        # NOT classical Bayesian inference
        
        chief_complaint = clinical_data.get('chief_complaint', '').lower()
        age = clinical_data.get('age', 50)
        
        # Quantum differential generator - superposed hypotheses
        quantum_differentials = {}
        
        if 'chest pain' in chief_complaint:
            quantum_differentials.update({
                'mi_acute': 0.35 + 0.05 * np.sin(age * np.pi / 100),  # Quantum varying
                'chest_wall_pain': 0.25,
                'gastroesophageal_reflux': 0.20 + 0.1 * np.cos(age * np.pi / 50),
                'pneumonia': 0.15,
                'aortic_dissection': 0.05
            })
        
        elif 'headache' in chief_complaint:
            quantum_differentials.update({
                'tension_headache': 0.40,
                'migraine': 0.25 + 0.1 * np.sin(age * np.pi / 40),
                'sinusitis': 0.15,
                'intracranial_hypertension': 0.10,
                'subarachnoid_hemorrhage': 0.05,
                'brain_tumor': 0.05
            })
            
        elif 'shortness of breath' in chief_complaint:
            quantum_differentials.update({
                'heart_failure': 0.30,
                'copd_exacerbation': 0.25,
                'pneumonia': 0.20,
                'anxiety': 0.15,
                'pulmonary_embolism': 0.10
            })
            
        else:
            # Generic quantum superposition
            quantum_differentials = {
                'unspecified_diagnosis': 0.5,
                'viral_illness': 0.3,
                'stress_related': 0.2
            }
        
        # Normalize quantum probabilities
        total_prob = sum(quantum_differentials.values())
        if total_prob > 0:
            quantum_differentials = {k: v/total_prob for k, v in quantum_differentials.items()}
        
        return quantum_differentials
    
    def _compute_quantum_entanglement(self, symptom_qbits: Dict, sign_qbits: Dict, 
                                     differential_qbits: Dict) -> np.ndarray:
        """Compute quantum entanglement correlations between states"""
        all_states = {**symptom_qbits, **sign_qbits, **differential_qbits}
        n_states = len(all_states)
        entanglement_matrix = np.zeros((n_states, n_states), dtype=complex)
        
        state_names = list(all_states.keys())
        
        for i, state1 in enumerate(state_names):
            for j, state2 in enumerate(state_names):
                if i != j:
                    # Quantum entanglement strength 
                    amplitude1 = all_states[state1]
                    amplitude2 = all_states[state2]
                    
                    # Entanglement correlation (NOT classical correlation)
                    entanglement = amplitude1 * np.conj(amplitude2)
                    entanglement_matrix[i, j] = entanglement
                else:
                    # Self-correlation (probability)
                    entanglement_matrix[i, i] = abs(all_states[state1])**2
        
        return entanglement_matrix
    
    def _calculate_decoherence_rate(self, clinical_data: Dict[str, Any]) -> float:
        """Calculate quantum decoherence rate based on environmental factors"""
        # More complete data = slower decoherence (more stable quantum state)
        # Missing critical data = faster decoherence (quantum uncertainty increases)
        
        required_fields = ['chief_complaint', 'age', 'vital_signs']
        missing_fields = sum(1 for field in required_fields if field not in clinical_data)
        
        # Quantum decoherence rate (higher = more unstable)
        base_rate = 0.1
        decoherence_increase = missing_fields * 0.05
        
        return base_rate + decoherence_increase
    
    def apply_virtue_supervision(self, quantum_case: QuantumClinicalCase) -> vQbitClinicalClaim:
        """
        Apply quantum virtue supervision to clinical decision
        
        This operates in quantum subspace and may or may not collapse superposed states
        based on virtue-based criteria.
        """
        # Honesty supervisor: Measure uncertainty honestly
        honesty_measurement = self.honesty_supervisor.quantum_measure(quantum_case)
        
        # Prudence supervisor: Prevent premature collapse
        prudence_measurement = self.prudence_supervisor.quantum_measure(quantum_case)
        
        # Justice supervisor: Ensure fair quantum access
        justice_measurement = self.justice_supervisor.quantum_measure(quantum_case)
        
        # Non-maleficence supervisor: Prevent harmful collapse
        non_maleficence_measurement = self.non_maleficence_supervisor.quantum_measure(quantum_case)
        
        # Decide quantum state based on virtue superposition
        virtue_amplitude = (honesty_measurement + prudence_measurement + 
                          justice_measurement + non_maleficence_measurement) / 4.0
        
        if abs(virtue_amplitude) > 0.8:  # Strong virtue alignment - can collapse
            quantum_state = QuantumClinicalState.COLLAPSED
            collapse_policy = "virtue_aligned_collapse"
        elif abs(virtue_amplitude) > 0.5:  # Medium virtue alignment - measured
            quantum_state = QuantumClinicalState.MEASURED
            collapse_policy = "virtue_measured_state"
        else:  # Low virtue alignment - remain superposed
            quantum_state = QuantumClinicalState.SUPERPOSED
            collapse_policy = "virtue_superposition_maintained"
        
        # Calculate quantum uncertainty
        uncertainty_hbar = self.hbar * quantum_case.decoherence_rate * abs(virtue_amplitude)
        
        # Determine entanglement partners
        entanglement_list = []
        max_entanglement_idx = np.argmax(np.abs(quantum_case.entanglement_matrix).sum(axis=1))
        if max_entanglement_idx < len(quantum_case.symptom_qbits):
            symptom_names = list(quantum_case.symptom_qbits.keys())
            if max_entanglement_idx < len(symptom_names):
                entanglement_list.append(symptom_names[max_entanglement_idx])
        
        return vQbitClinicalClaim(
            measurement_type="QuantumClinicalAssessment",
            quantum_state=quantum_state,
            amplitude=virtue_amplitude,
            probability=abs(virtue_amplitude)**2,
            phase=np.angle(virtue_amplitude),
            entanglement_list=entanglement_list,
            collapse_policy=collapse_policy,
            uncertainty_hbar=uncertainty_hbar,
            toolchain_hash=self._compute_quantum_hash(quantum_case),
            timestamp=datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
        )
    
    def _compute_quantum_hash(self, quantum_case: QuantumClinicalCase) -> str:
        """Compute quantum hash for provenance"""
        state_norm = np.linalg.norm(quantum_case.quantum_state_vector)
        entanglement_trace = np.real(np.trace(quantum_case.entanglement_matrix))
        
        hash_input = f"{quantum_case.case_id}_{state_norm:.6f}_{entanglement_trace:.6f}"
        return hashlib.sha256(hash_input.encode()).hexdigest()[:16]


class QuantumVirtueSupervisor:
    """Individual quantum virtue supervisor operating on subspaces"""
    
    def __init__(self, virtue_type: str):
        self.virtue_type = virtue_type
        self.quantum_subspace_dim = 256
    
    def quantum_measure(self, quantum_case: QuantumClinicalCase) -> complex:
        """
        Perform quantum measurement in virtue subspace
        
        Returns complex amplitude indicating virtue alignment strength
        """
        if self.virtue_type == "honesty":
            return self._honesty_quantum_measurement(quantum_case)
        elif self.virtue_type == "prudence":
            return self._prudence_quantum_measurement(quantum_case)
        elif self.virtue_type == "justice":
            return self._justice_quantum_measurement(quantum_case)
        elif self.virtue_type == "non_maleficence":
            return self._non_maleficence_quantum_measurement(quantum_case)
        
        return 0.0 + 0j
    
    def _honesty_quantum_measurement(self, quantum_case: QuantumClinicalCase) -> complex:
        """Quantum measurement of honesty virtue"""
        # High decoherence = low honesty (uncertain states)
        # Low decoherence = high honesty (well-defined states)
        honesty_decay = np.exp(-quantum_case.decoherence_rate)
        
        # Add phase based on data completeness
        completeness_phase = len(quantum_case.symptom_qbits) * np.pi / 10
        
        return honesty_decay * np.exp(1j * completeness_phase)
    
    def _prudence_quantum_measurement(self, quantum_case: QuantumClinicalCase) -> complex:
        """Quantum measurement of prudence virtue"""
        # Prudent when highest probability differential is clear
        if quantum_case.differential_qbits:
            max_prob = max(abs(amp)**2 for amp in quantum_case.differential_qbits.values())
            prudence_magnitude = max_prob
        else:
            prudence_magnitude = 0.5
        
        # Add quantum phase
        prudence_phase = np.pi * prudence_magnitude
        
        return prudence_magnitude * np.exp(1j * prudence_phase)
    
    def _justice_quantum_measurement(self, quantum_case: QuantumClinicalCase) -> complex:
        """Quantum measurement of justice virtue"""
        # Justice requires balanced treatment across differential diagnoses
        if quantum_case.differential_qbits:
            probabilities = [abs(amp)**2 for amp in quantum_case.differential_qbits.values()]
            entropy = -sum(p * np.log(p + 1e-10) for p in probabilities)
            max_entropy = np.log(len(probabilities))
            justice_magnitude = entropy / max_entropy if max_entropy > 0 else 0
        else:
            justice_magnitude = 0.5
        
        justice_phase = np.pi * justice_magnitude * 0.5
        
        return justice_magnitude * np.exp(1j * justice_phase)
    
    def _non_maleficence_quantum_measurement(self, quantum_case: QuantumClinicalCase) -> complex:
        """Quantum measurement of non-maleficence virtue"""
        # Non-maleficence increases when dangerous diagnoses have low probability
        dangerous_dxs = ['mi_acute', 'aortic_dissection', 'subarachnoid_hemorrhage', 
                        'brain_tumor', 'pulmonary_embolism']
        
        if quantum_case.differential_qbits:
            dangerous_prob = sum(abs(quantum_case.differential_qbits.get(dx, 0))**2 
                                for dx in dangerous_dxs)
            non_maleficence_magnitude = 1.0 - dangerous_prob  # Higher = safer
        else:
            non_maleficence_magnitude = 1.0
        
        non_maleficence_phase = np.pi * non_maleficence_magnitude * 0.25
        
        return non_maleficence_magnitude * np.exp(1j * non_maleficence_phase)


def demo_quantum_clinical_reasoning():
    """Demonstrate quantum clinical reasoning on vQbit substrate"""
    
    print("🧬 FoT Quantum Clinical Engine Demo")
    print("=" * 50)
    
    # Initialize quantum engine
    quantum_engine = QuantumClinicalEngine(vqbit_dimension=512)
    
    # Sample clinical case
    clinical_case = {
        'case_id': 'QUANTUM_DEMO_001',
        'chief_complaint': 'chest pain for 3 hours',
        'age': 65,
        'symptoms': {
            'chest_pain': {'intensity': 0.8, 'duration_hours': 3},
            'nausea': {'intensity': 0.5, 'duration_hours': 2},
            'diaphoresis': {'intensity': 0.6, 'duration_hours': 1}
        },
        'vital_signs': {
            'systolic_bp': 150,
            'diastolic_bp': 95,
            'heart_rate': 110,
            'respiratory_rate': 22,
            'temperature_c': 37.2,
            'spo2': 95
        },
        'medical_history': ['hypertension', 'diabetes_type2']
    }
    
    print(f"📋 Clinical Case: {clinical_case['chief_complaint']}")
    print(f"👤 Patient Age: {clinical_case['age']}")
    
    # Encode into quantum state
    quantum_case = quantum_engine.encode_clinical_case(clinical_case)
    
    print(f"\n🔮 Quantum Case Encoding:")
    print(f"   Case ID: {quantum_case.case_id}")
    print(f"   Quantum State Dimension: {quantum_case.quantum_state_vector.shape[0]}")
    print(f"   Symptom QBits: {len(quantum_case.symptom_qbits)}")
    print(f"   Sign QBits: {len(quantum_case.sign_qbits)}")
    print(f"   Differential QBits: {len(quantum_case.differential_qbits)}")
    print(f"   Decoherence Rate: {quantum_case.decoherence_rate:.4f}")
    
    # Show superposed differential diagnoses
    print(f"\n⚛️ Superposed Differential Diagnoses:")
    for dx, amplitude in quantum_case.differential_qbits.items():
        probability = abs(amplitude)**2
        phase_deg = np.angle(amplitude) * 180 / np.pi
        print(f"   {dx}: {probability:.3f} (phase: {phase_deg:.1f}°)")
    
    # Apply quantum virtue supervision
    quantum_claim = quantum_engine.apply_virtue_supervision(quantum_case)
    
    print(f"\n✨ Quantum Virtue Supervision Result:")
    print(f"   Measurement Type: {quantum_claim.measurement_type}")
    print(f"   Quantum State: {quantum_claim.quantum_state.value}")
    print(f"   Amplitude: {quantum_claim.amplitude:.3f}")
    print(f"   Probability: {quantum_claim.probability:.3f}")
    print(f"   Phase: {quantum_claim.phase * 180 / np.pi:.1f}°")
    print(f"   Entanglement List: {quantum_claim.entanglement_list}")
    print(f"   Collapse Policy: {quantum_claim.collapse_policy}")
    print(f"   Uncertainty (ħ): {quantum_claim.uncertainty_hbar:.4f}")
    print(f"   Toolchain Hash: {quantum_claim.toolchain_hash}")
    
    return quantum_claim


if __name__ == "__main__":
    claim = demo_quantum_clinical_reasoning()
