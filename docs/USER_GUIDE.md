# 👥 FoTClinician User Guide

## 🎯 Complete User Manual for Quantum Clinical Decision Support

---

## 📋 **Table of Contents**

1. [🚀 Getting Started](#-getting-started)
2. [🔮 Quantum Diagnosis Interface](#-quantum-diagnosis-interface)  
3. [📊 Data Readiness Validation](#-data-readiness-validation)
4. [🛡️ Virtue Supervision Panel](#️-virtue-supervision-panel)
5. [📋 Medical Coding](#-medical-coding)
6. [⚙️ Configuration](#️-configuration)
7. [❓ Troubleshooting](#-troubleshooting)

---

## 🚀 **Getting Started**

### Prerequisites

Before using FoTClinician, ensure you have:

- ✅ **Python 3.8+** installed
- ✅ **Medical knowledge** or clinical background
- ✅ **Web browser** for Streamlit interface
- ✅ **Valid clinical data** for analysis

### Installation

```bash
# Clone the repository
git clone https://github.com/FortressAI/FoTClinician.git
cd FoTClinician

# Install dependencies
pip install -r requirements.txt

# Launch the application
streamlit run streamlit_app.py
```

### First Launch

1. **🌐 Open your browser** to `http://localhost:8501`
2. **🎯 Navigate** to "Quantum Diagnosis Demo"
3. **🚀 Click** "MI Case Template" for demonstration
4. **⚡ Execute** "Quantum Clinical Analysis"

---

## 🔮 **Quantum Diagnosis Interface**

### Overview

The quantum diagnosis interface demonstrates **vQbit quantum substrate** clinical reasoning with **virtue-based decision making**.

### Case Input Methods

#### 🚀 **Quick Case Templates**

**Available Templates:**
- 🏥 **MI Case**: Acute myocardial infarction scenario
- 🧸 **Pediatric Case**: Febrile infant sepsis example

**Usage:**
1. Click desired template button
2. Review auto-populated fields  
3. Modify as needed
4. Execute analysis

#### ✍️ **Manual Case Entry**

**Required Fields:**
```json
{
  "case_id": "Unique identifier",
  "chief_complaint": "Primary medical concern",
  "age": "45",
  "gender": "male|female|unknown"
}
```

**Optional Fields:**
```json
{
  "medical_history": ["condition1", "condition2"],
  "vital_signs": {
    "systolic_bp": 120,
    "diastolic_bp": 80,
    "heart_rate": 75,
    "temperature_c": 37.0
  },
  "symptoms": {
    "chest_pain": {"intensity": 0.8},
    "shortness_breath": {"intensity": 0.6}
  }
}
```

### Understanding Quantum Results

#### ⚛️ **Quantum Superposed Diagnoses**

**Probability Distribution Chart:**
- **Red bars**: High-probability diagnoses
- **Blue bars**: Lower-probability alternatives  
- **Height**: Quantum probability amplitude

**Diagnosis List:**
- **Quantum Probability**: Probability of each hypothesis
- **Quantum Amplitude**: Complex amplitude magnitude
- **Evidence Strength**: Diagnostic confidence

#### 🛡️ **Virtue Supervision Results**

**Key Metrics:**
- **Quantum State**: `superspeed` → `measured` → `collapsed`
- **Uncertainty**: Higher values indicate more uncertain cases
- **Collapse Policy**: Virtue constraint applied
- **QBits Count**: Number of quantum states considered

---

## 📊 **Data Readiness Validation**

### Purpose

Ensures clinical cases have sufficient data for reliable quantum diagnosis.

### Validation Tracks

#### 🔍 **Medication Safety Track**
*Requires: medication list, allergies, renal/hepatic function*

**Validation Criteria:**
- ✅ Complete medication inventory
- ✅ Allergy list documented  
- ✅ Renal function (eGFR/creatinine <12 months)
- ✅ Hepatic panel (AST, ALT, bilirubin <12 months)
- ✅ Pregnancy status (females 12-55 years)

#### 🚨 **Triage Assessment Track** 
*Requires: chief complaint + 5 vital signs*

**Validation Criteria:**
- ✅ Chief complaint documented
- ✅ Blood pressure (systolic/diastolic)
- ✅ Heart rate
- ✅ Respiratory rate
- ✅ Temperature
- ✅ Oxygen saturation (SpO2)

#### 🔬 **Next Diagnostic Step Track**
*Requires: complaint + clinical context*

**Validation Criteria:**
- ✅ Chief complaint documented
- ✅ Clinical context (HPI OR active problems OR labs/imaging)

### Validation Interface

**🔍 Comprehensive Case Validation:**
1. Input clinical case JSON
2. Click "Validate Data Readiness"
3. Review track-by-track results
4. Address identified gaps

**⚡ Quick Validation:**
1. Check required data components
2. View completeness percentage
3. Get readiness status

---

## 🛡️ **Virtue Supervision Panel**

### Cardinal Virtues

#### 🎯 **Honesty Virtue**
**Function**: Surface uncertainty genuinely
**Quantum Application**: Prevents false confidence in diagnoses
**Monitoring**: Uncertainty transparency metrics

#### 🤔 **Prudence Virtue**  
**Function**: Default to safest medical options
**Quantum Application**: Conservative diagnostic collapses
**Monitoring**: Safety-first decision protocols

#### ⚖️ **Justice Virtue**
**Function**: Fair resource allocation
**Quantum Application**: Bias prevention algorithms
**Monitoring**: Fairness distribution metrics

#### 🚫 **Non-maleficence Virtue**
**Function**: Prevent harm to patients
**Quantum Application**: Harm-blocking quantum gates
**Monitoring**: Safety constraint enforcement

### Monitoring Dashboard

**Virtue Scores Radar Chart:**
- Shows real-time virtue compliance
- Alerts on virtue conflicts
- Tracks ethical decision patterns

**Recent Conflicts Log:**
- Timestamp of virtue conflicts
- Conflict type identification
- Resolution approach taken
- Quantum system impact assessment

---

## 📋 **Medical Coding**

### CPT (Current Procedural Terminology)

#### 🩺 **Automatic Code Generation**

**Process:**
1. Select primary diagnoses
2. System generates relevant CPT codes
3. Review recommendations
4. Export for billing

**Code Information:**
- **Code**: 5-digit CPT identifier
- **Description**: Procedural details
- **RVU**: Relative Value Units
- **Medicare Fee**: Estimated reimbursement

#### 🩺 **Supported Procedures**

**Cardiac Procedures:**
- EKG Interpretation (93010): $15.00
- Cardiac Catheterization (93458): $52.00

**Laboratory Tests:**
- Blood Glucose (82948): $23.00
- Basic Metabolic Panel (80053): $16.00

### ICD-10 (Diagnostic Codes)

#### 📊 **Diagnosis Mapping**

**Common Diagnoses:**
- Acute MI: `I21.9`
- Diabetic Ketoacidosis: `E10.10`
- Acute Appendicitis: `K35.80`
- Hypertensive Crisis: `I16.9`
- Sepsis: `A41.9`

#### 📊 **Usage Process**
1. Confirm quantum diagnosis
2. System maps to ICD-10 code
3. Review description accuracy
4. Validate billing appropriateness

### DRG (Diagnosis Related Groups)

#### 💰 **Payment Analysis**

**Complexity Levels:**
- **Low**: Simple cases, lower reimbursement
- **Medium**: Standard complexity
- **High**: Complex cases requiring expertise  
- **Critical**: Intensive care scenarios

**Service Types:**
- **Medical**: Non-surgical management
- **Surgical**: Operative procedures
- **Emergency**: Urgent/emergent cases
- **Critical Care**: ICU-level management

**DRG Information:**
- **Code**: 3-digit DRG identifier
- **Base Payment**: Medicare reimbursement amount
- **Expected LOS**: Length of stay estimate

---

## ⚙️ **Configuration**

### Quantum Engine Settings

**File**: `core/clinical/quantum_clinical_engine.py`

```python
# Quantum substrate configuration
vqbit_dimension = 512          # Quantum state space size
quantum_decoherence_rate = 0.1  # Natural uncertainty parameter
virtue_supervisor_enabled = True # Ethical constraint activation
```

### Validation Settings

**File**: `core/clinical/data_readiness_checker.py`

```python
# Medical validation thresholds
passing_threshold = 0.85        # 85% required for medical grade
safety_protocol_level = "maximum"  # Maximum safety enforcement
usmle_validation_mode = "strict"   # Rigorous medical standards
```

### Interface Settings

**File**: `streamlit_app.py`

```python
# UI configuration
page_title = "FoTClinician - Quantum Medical AI"
page_icon = "🧠⚛️"
layout = "wide"
initial_sidebar_state = "expanded"
```

---

## ❓ **Troubleshooting**

### Common Issues

#### 🔧 **"Module not found" errors**
```bash
# Solution: Install dependencies
pip install -r requirements.txt
```

#### 🌐 **Streamlit not launching**
```bash
# Solution: Check port availability
streamlit run streamlit_app.py --server.port 8502
```

#### ⚛️ **Quantum analysis fails**
```bash
# Solution: Validate input data format
python -c "import json; json.loads('your_case_data')"
```

#### 🧪 **USMLE tests failing**
```bash
# Solution: Run specific validation
python tests/test_usmle_board_certification.py
```

### Error Messages

#### ❌ **"Quantum analysis error"**
- **Cause**: Invalid input data format
- **Solution**: Verify JSON structure, required fields present

#### ❌ **"Data readiness validation failed"**  
- **Cause**: Insufficient clinical data
- **Solution**: Add missing required fields, increase completeness

#### ❌ **"Virtue supervisor conflict"**
- **Cause**: Ethical constraint violation
- **Solution**: Review case safety, modify inputs if appropriate

### Support Resources

#### 📞 **Getting Help**
- 🐛 **Report Issues**: [GitHub Issues](https://github.com/FortressAI/FoTClinician/issues)
- 💬 **Community Forum**: [GitHub Discussions](https://github.com/FortressAI/FoTClinician/discussions)
- 📧 **Direct Contact**: fortress.ai@clinical.support

#### 📖 **Additional Documentation**
- 🎓 **[USMLE Certification Guide](USMLE_CERTIFICATION.md)**: Medical exam validation
- ⚛️ **[Quantum Reference](QUANTUM_REFERENCE.md)**: Scientific methodology
- 🛡️ **[Safety Protocols](SAFETY_GUIDE.md)**: Non-maleficence procedures

---

## ⚖️ **Legal & Safety Information**

### Medical Disclaimer

⚠️ **IMPORTANT**: FoTClinician is a **decision-support tool only**. It cannot replace professional medical judgment or qualified healthcare providers.

### Usage Guidelines

1. 🔄 **Always validate** quantum diagnoses with clinical expertise
2. ⚠️ **Emergency situations**: Contact emergency services immediately  
3. 🔒 **Data privacy**: Ensure patient data protection compliance
4. 📋 **Documentation**: Maintain audit trails of AI recommendations

### Ethical Considerations

- ⚖️ **Fair Use**: Avoid algorithmic bias in patient selection
- 🔒 **Privacy**: Protect patient identifying information
- 🤝 **Collaboration**: Work with qualified medical professionals
- 📊 **Transparency**: Document AI-assisted decision making

---

*📚 This guide covers version 1.0 of FoTClinician Quantum Clinical Decision Support System*

**© 2024 Fortress AI - Quantum Medical Intelligence**
