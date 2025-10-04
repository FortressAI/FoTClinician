# 🔌 **API Reference - Core Clinical Engine**

**Technical Documentation for Quantum Clinical AI Integration**

---

## 🎯 **API Overview**

FoTClinician provides programmatic access to quantum clinical intelligence through a comprehensive REST API designed for EMR integration, clinical decision support automation, and healthcare workflow optimization.

---

## 🚀 **Authentication**

### **🔐 API Key Management**

**Getting Started:**
```bash
# Request API access
curl -X POST https://api.fotclinician.com/v1/auth/request \
  -H "Content-Type: application/json" \
  -d '{
    "organization": "Your Hospital Name",
    "contact_email": "admin@your-hospital.com",
    "intended_use": "EMR Integration"
  }'
```

**Authentication Header:**
```http
GET /v1/clinical/analyze
Authorization: Bearer YOUR_API_KEY
Content-Type: application/json
```

### **🔑 Permissions & Scopes**

**Available Scopes:**
- `clinical:analyze` - Run quantum diagnostic analysis
- `coding:generate` - Generate medical codes
- `validation:check` - Data readiness validation
- `virtue:monitor` - Ethical constraint monitoring

---

## 🔮 **Quantum Clinical Analysis API**

### **🔬 Core Analysis Endpoint**

**Endpoint:** `POST /v1/clinical/analyze`

**Request Schema:**
```json
{
  "case_id": "string",
  "patient_data": {
    "chief_complaint": "string",
    "age": "integer",
    "gender": "string",
    "vital_signs": {
      "systolic_bp": "number",
      "diastolic_bp": "number", 
      "heart_rate": "number",
      "respiratory_rate": "number",
      "temperature_c": "number",
      "spo2": "number"
    },
    "symptoms": "object",
    "medical_history": "array",
    "medications": "array",
    "allergies": "array",
    "laboratory": "object"
  },
  "specialty_focus": "string",
  "urgency_level": "string",
  "quantum_params": {
    "vqbit_dimension": 512,
    "virtue_thresholds": {
      "honesty": 0.25,
      "prudence": 0.30,
      "justice": 0.20,
      "non_maleficence": 0.35
    }
  }
}
```

**Response Schema:**
```json
{
  "analysis_id": "string",
  "timestamp": "ISO8601",
  "quantum_case": {
    "case_id": "string",
    "quantum_state_dimension": 512,
    "diagnostic_qbits": {
      "diagnosis_name": {
        "amplitude": "number",
        "phase": "number",
        "probability": "number"
      }
    },
    "decay_rate": "number",
    "entanglement_list": ["string"],
    "decoherence_time": "number"
  },
  "clinical_recommendations": {
    "primary_diagnosis": {
      "name": "string",
      "confidence": "number",
      "evidence_support": ["string"],
      "next_steps": ["string"]
    },
    "differential_diagnoses": [
      {
        "name": "string",
        "probability": "number",
        "clinical_indicators": ["string"]
      }
    ]
  },
  "virtue_supervision": {
    "measurement_type": "string",
    "quantum_state": "measured|superposed",
    "amplitude": "complex_number",
    "phase": "number",
    "probability": "number",
    "uncertainty": "number",
    "collapsed_disorder": "string",
    "virtue_constraints": {
      "honesty": "applied|violated|pending",
      "prudence": "applied|violated|pending",
      "justice": "applied|violated|pending", 
      "non_maleficence": "applied|violated|pending"
    },
    "toolchain_hash": "string"
  },
  "usmle_validation": {
    "step_1_compliance": "number",
    "step_2_ck_compliance": "number",
    "step_3_compliance": "number",
    "overall_certification": "certified|pending|failed"
  },
  "metadata": {
    "processing_time_ms": "integer",
    "quantum_operations": "integer",
    "confidence_threshold": "number",
    "specialty_expertise": "string"
  }
}
```

### **🎯 Example Usage**

**Python Integration:**
```python
import requests
import json

def analyze_clinical_case(patient_data, api_key):
    """
    Analyze clinical case using FoTClinician quantum engine
    """
    url = "https://api.fotclinician.com/v1/clinical/analyze"
    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json"
    }
    
    payload = {
        "case_id": f"quantum_analysis_{timestamp()}",
        "patient_data": patient_data,
        "specialty_focus": "emergency_medicine",
        "urgency_level": "urgent",
        "quantum_params": {
            "vqbit_dimension": 512,
            "virtue_thresholds": {
                "honesty": 0.25,
                "prudence": 0.30,
                "non_maleficence": 0.35
            }
        }
    }
    
    response = requests.post(url, headers=headers, json=payload)
    
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"API Error: {response.status_code} - {response.text}")

# Example cardiac emergency case
cardiac_case = {
    "chief_complaint": "crushing chest pain radiating to left arm",
    "age": 55,
    "gender": "male",
    "vital_signs": {
        "systolic_bp": 160,
        "diastolic_bp": 110,
        "heart_rate": 110,
        "temperature_c": 37.0,
        "spo2": 94
    },
    "medical_history": ["hypertension", "smoking"],
    "current_medications": ["lisinopril"],
    "ecg_findings": {
        "st_elevation_anterior": True,
        "q_waves": False
    },
    "laboratory": {
        "troponin_i": 4.2,
        "ck_mb": 95
    }
}

result = analyze_clinical_case(cardiac_case, "YOUR_API_KEY")
```

---

## 📋 **Medical Coding API**

### **🏥 Automatic Code Generation**

**Endpoint:** `POST /v1/coding/generate`

**Request Schema:**
```json
{
  "analysis_result": "object",
  "specialty": "string",
  "encounter_type": "string",
  "coding_preferences": {
    "include_cpt": true,
    "include_icd10": true,
    "include_drg": true,
    "validate_codes": true
  }
}
```

**Response Schema:**
```json
{
  "codes": {
    "primary_diagnosis": {
      "icd10": "string",
      "description": "string",
      "confidence": "number"
    },
    "procedures": [{
      "cpt": "string",
      "description": "string",
      "rvu": "number",
      "modifiers": ["string"]
    }],
    "complications": [{
      "icd10": "string", 
      "description": "string",
      "severity": "string"
    }],
    "secondary_diagnoses": [{
      "icd10": "string",
      "description": "string",
      "chronicity": "acute|chronic|subacute"
    }]
  },
  "reimbursement": {
    "primary_drg": "string",
    "mdc": "string", 
    "case_mix_index": "number",
    "estimated_payment": "number",
    "complexity_level": "string"
  },
  "validation": {
    "code_accuracy": "number",
    "compliance_check": "passed|failed|warning",
    "documentation_gaps": ["string"]
  }
}
```

---

## 📊 **Data Readiness Validation API**

### **🕘 Validation Check**

**Endpoint:** `POST /v1/validation/check`

**Request Schema:**
```json
{
  "clinical_data": "object",
  "validation_tracks": ["medication_safety", "triage_assessment", "diagnostic_planning"],
  "required_thresholds": {
    "medication_safety": 1.0,
    "triage_assessment": 1.0,
    "diagnostic_planning": 1.0
  }
}
```

**Response Schema:**
```json
{
  "overall_readiness": "ready|not_ready|near_miss",
  "track_results": [{
    "track_name": "string",
    "result": "ready|not_ready|near_miss",
    "score": "number",
    "gaps": [{
      "field": "string",
      "reason": "string", 
      "severity": "required|highly_recommended|optional",
      "example_value": "any"
    }],
    "recommendations": ["string"]
  }],
  "collapsed_ready": "boolean",
  "validation_metadata": {
    "check_timestamp": "ISO8601",
    "data_provenance": "string",
    "contract_version": "string"
  }
}
```

---

## 🛡️ **Virtue Supervision API**

### **⚡ Ethical Constraint Monitoring**

**Endpoint:** `POST /v1/virtue/supervise`

**Request Schema:**
```json
{
  "quantum_state": "object",
  "clinical_context": {
    "urgency": "routine|urgent|emergency|critical",
    "patient_age": "number",
    "coercion_present": "boolean",
    "resource_constraints": "object"
  },
  "virtue_settings": {
    "honesty_threshold": "number",
    "prudence_threshold": "number", 
    "justice_threshold": "number",
    "non_maleficence_threshold": "number"
  }
}
```

**Response Schema:**
```json
{
  "supervision_result": {
    "honesty": {
      "applied": "boolean",
      "uncertainty_surfaced": "boolean",
      "confidence_reporting": "comprehensive|limited"
    },
    "prudence": {
      "applied": "boolean",
      "safety_focused": "boolean", 
      "conservative_default": "boolean"
    },
    "justice": {
      "applied": "boolean",
      "bias_detected": "boolean",
      "fair_treatment": "boolean"
    },
    "non_maleficence": {
      "applied": "boolean",
      "harms_prevented": ["string"],
      "contraindications_flaggeed": ["string"]
    }
  },
  "ethical_debt": "number",
  "monitoring_alerts": ["string"],
  "compliance_status": "compliant|warning|violation"
}
```

---

## ⚡ **Rate Limits & Performance**

### **📈 API Limits**

**Tier Structure:**
- **Free Tier**: 100 requests/day
- **Professional**: 10,000 requests/day  
- **Enterprise**: Unlimited with SLA

**Rate Limiting:**
```
X-RateLimit-Limit: 10000
X-RateLimit-Remaining: 9999
X-RateLimit-Reset: 1647129600
```

### **🚀 Performance SLA**

**Response Times:**
- P95 Response Time: <2.0s
- P99 Response Time: <5.0s
- Availability: 99.9%

**Throughput:**
- Max Concurrent Requests: 1000
- Requests per Second: 500
- Burst Capacity: 1000 req/s for 10s

---

## 🔍 **Error Handling**

### **❌ Standard Error Responses**

**HTTP Status Codes:**
- `200` - Success
- `400` - Bad Request (invalid data)
- `401` - Unauthorized (invalid API key)
- `402` - Payment Required (quota exceeded)
- `403` - Forbidden (insufficient permissions)
- `422` - Validation Error (invalid parameters)
- `429` - Rate Limited (too many requests)
- `500` - Internal Server Error
- `503` - Service Unavailable (maintenance)

**Error Response Schema:**
```json
{
  "error": {
    "code": "string",
    "message": "string",
    "details": {
      "field": "string",
      "violation": "string",
      "suggested_fix": "string"
    },
    "timestamp": "ISO8601",
    "request_id": "string"
  }
}
```

### **🔧 Common Error Codes**

**Validation Errors:**
- `INVALID_PATIENT_DATA` - Malformed clinical data
- `MISSING_REQUIRED_FIELD` - Required patient information absent
- `INVALID_VITAL_SIGNS` - Unrealistic vital sign values
- `INCOMPATIBLE_STATUS` - Conflicting clinical status indicators

**Processing Errors:**
- `QUANTUM_PROCESSING_FAILED` - Internal quantum computation error
- `VIRTUE_CONSTRAINT_VIOLATION` - Ethical constraint enforcement failure
- `USMLE_VALIDATION_ERROR` - Clinical validation algorithm error
- `CODING_GENERATION_FAILED` - Medical code generation error

---

## 🔄 **Webhook Integration**

### **📡 Event Notifications**

**Webhook Configuration:**
```bash
curl -X POST https://api.fotclinician.com/v1/webhooks \
  -H "Authorization: Bearer YOUR_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "url": "https://your-system.com/webhooks/fotclinician",
    "events": ["analysis.completed", "virtue.alert", "validation.failed"],
    "secret": "your_webhook_secret"
  }'
```

**Webhook Payload:**
```json
{
  "event_type": "analysis.completed",
  "timestamp": "ISO8601",
  "data": {
    "analysis_id": "string",
    "case_id": "string", 
    "status": "completed|failed|warning",
    "results": "object"
  },
  "signature": "HMAC_SHA256_SIGNATURE"
}
```

---

## 📚 **SDKs & Libraries**

### **🐍 Python SDK**

**Installation:**
```bash
pip install fotclinician
```

**Usage:**
```python
from fotclinician import ClinicalAI

client = ClinicalAI(api_key="YOUR_API_KEY")

# Analyze clinical case
result = client.analyze(
    patient_data=cardiac_case,
    specialty="emergency_medicine"
)

# Generate medical codes
codes = client.generate_codes(
    analysis_result=result,
    specialty="cardiology"
)

# Validate data readiness
validation = client.validate_readiness(
    clinical_data=case_data
)
```

### **🟨 JavaScript SDK**

**Installation:**
```bash
npm install @fortress-ai/fotclinician
```

**Usage:**
```javascript
import { ClinicalAI } from '@fortress-ai/fotclinician';

const client = new ClinicalAI('YOUR_API_KEY');

// Analyze clinical case
const result = await client.analyze({
  patientData: cardiacCase,
  specialty: 'emergency_medicine'
});

// Generate codes
const codes = await client.generateCodes({
  analysisResult: result,
  specialty: 'cardiology'
});
```

---

## 📖 **Testing & Development**

### **🧪 Sandbox Environment**

**Test API Key:** `sb_[32_character_preview]`

**Test Endpoints:**
- Analysis: `https://api-sandbox.fotclinician.com/v1/clinical/analyze`
- Coding: `https://api-sandbox.fotclinician.com/v1/coding/generate`
- Validation: `https://api-sandbox.fotclinician.com/v1/validation/check`

### **🔨 Development Tools**

**API Explorer:**
- URL: `https://developers.fotclinician.com/explorer`
- Features: Interactive API testing, documentation browser
- Authentication: Test token integration

**Mock Server:**
```bash
# Start local mock server
npm install -g @fortress-ai/mock-server
mock-server start --port 3001

# Test endpoints locally
curl http://localhost:3001/v1/clinical/analyze \
  -X POST \
  -H "Content-Type: application/json" \
  -d @test_data.json
```

---

## 📞 **Support & Resources**

### **🆘 Developer Support**

**Documentation:**
- [Interactive API Explorer](https://api.fotclinician.com/explorer)
- [Python SDK Docs](https://docs.fotclinician.com/python)
- [JavaScript SDK Docs](https://docs.fotclinician.com/javascript)
- [Webhook Guide](https://docs.fotclinician.com/webhooks)

**Community:**
- [GitHub Repository](https://github.com/FortressAI/FoTClinician)
- [Discord Community](https://discord.gg/fotclinician)
- [Stack Overflow](https://stackoverflow.com/search?q=fotclinician)

**Technical Support:**
- Email: developers@fortress.ai
- Discord: Real-time developer chat
- GitHub Issues: Bug reports and feature requests

---

**🔌 Ready to integrate quantum clinical intelligence into your applications?**

**Start building with our comprehensive API today and revolutionize healthcare workflows with AI!**
