# 🖼️🎵 Medical Data Quality Validation

**Comprehensive Image and Audio Readiness Framework for Clinical Data**

---

## 🎯 **Overview**

FoTClinician includes a comprehensive medical data quality validation framework that ensures clinical images and audio meet professional medical standards before processing by the quantum clinical engine.

---

## 🖼️ **Image Readiness Validation**

### **Core Functionality**
The image readiness checker validates medical imaging data against Field of Truth clinical readiness requirements, checking for PHI compliance, quality metrics, and technical specifications.

### **Quality Metrics Validated**
- **Focus Score**: Image sharpness and clarity (threshold: >0.6)
- **Exposure Score**: Proper lighting and contrast (range: 0.3-0.9)
- **Signal-to-Noise Ratio**: Image quality (threshold: >20 dB)
- **Resolution**: Minimum dimensions (threshold: >512px shortest side)

### **Technical Requirements**
- **Modality**: Imaging type (fundus, derm, X-ray, etc.)
- **Body Site**: Anatomical location
- **Pixel Spacing**: Calibrated measurements OR scale reference
- **Device Model**: Acquisition equipment information
- **Acquisition Time**: Timestamp of image capture

### **PHI Compliance**
- **PHI Burn-in Flag**: Automatic detection of patient identifiers
- **Privacy Protection**: Blocks processing if PHI detected
- **De-identification**: Ensures HIPAA compliance

### **Test Cases**

#### **✅ Good Quality Example: Retina Fundus**
```json
{
  "modality": "fundus",
  "bodySite": "retina-left",
  "widthPx": 1024,
  "heightPx": 1024,
  "pixelSpacingMm": 0.012,
  "qualityMeasurements": [
    {"hasMetric": "fimg:Quality_FocusScore", "value": 0.81},
    {"hasMetric": "fimg:Quality_ExposureScore", "value": 0.60},
    {"hasMetric": "fimg:Quality_SNR_dB", "value": 26.5}
  ]
}
```
**Result**: ✅ READY - All quality metrics pass

#### **❌ Near-Miss Example: Dermatology**
```json
{
  "modality": "derm",
  "bodySite": "skin-forearm",
  "widthPx": 480,
  "heightPx": 360,
  "pixelSpacingMm": null,
  "qualityMeasurements": [
    {"hasMetric": "fimg:Quality_FocusScore", "value": 0.44},
    {"hasMetric": "fimg:Quality_ExposureScore", "value": 0.95},
    {"hasMetric": "fimg:Quality_SNR_dB", "value": 18.0}
  ]
}
```
**Result**: ❌ NOT READY - Poor focus, overexposed, low SNR, insufficient resolution

---

## 🎵 **Audio Readiness Validation**

### **Core Functionality**
The audio readiness checker validates medical audio data against clinical acquisition standards, ensuring proper calibration, quality metrics, and technical specifications.

### **Quality Metrics Validated**
- **Signal-to-Noise Ratio**: Audio clarity (threshold: >20 dB)
- **Noise Floor**: Background noise level (threshold: <-35 dBFS)
- **Artifact Score**: Signal distortion (threshold: <0.4)
- **Calibration Status**: Device calibration verification

### **Technical Requirements**
- **Body Site**: Anatomical location (heart, lung, etc.)
- **Sample Rate**: Audio frequency (threshold: >4000 Hz)
- **Bit Depth**: Audio resolution (threshold: >16 bits)
- **Channels**: Mono audio preferred (warning if stereo)
- **Duration**: Recording length (recommended: >10 seconds)

### **Clinical Standards**
- **Heart Sounds**: Cardiac auscultation quality
- **Lung Sounds**: Respiratory assessment quality
- **Environmental**: Quiet room conditions
- **Patient Position**: Proper positioning for acquisition

### **Test Cases**

#### **✅ Good Quality Example: Lung Audio**
```json
{
  "bodySite": "posterior-lower-left",
  "sampleRateHz": 8000,
  "bitDepth": 16,
  "channels": 1,
  "durationSec": 20.0,
  "calibrationPassed": 1,
  "qualityMeasurements": [
    {"hasMetric": "faud:Quality_SNR_dB", "value": 24.0},
    {"hasMetric": "faud:Quality_NoiseFloor_dBFS", "value": -42.0},
    {"hasMetric": "faud:Quality_ArtifactScore", "value": 0.18}
  ]
}
```
**Result**: ✅ READY - All quality metrics pass

#### **❌ Near-Miss Example: Heart Audio**
```json
{
  "bodySite": "apex",
  "sampleRateHz": 4000,
  "bitDepth": 16,
  "channels": 1,
  "durationSec": 8.0,
  "calibrationPassed": 0,
  "qualityMeasurements": [
    {"hasMetric": "faud:Quality_SNR_dB", "value": 12.0},
    {"hasMetric": "faud:Quality_NoiseFloor_dBFS", "value": -25.0},
    {"hasMetric": "faud:Quality_ArtifactScore", "value": 0.52}
  ]
}
```
**Result**: ❌ NOT READY - Low SNR, high noise floor, calibration failed, short duration

---

## 🧪 **Synthetic Test Data**

### **Generated Test Files**
The synthetic data generator creates realistic test files for immediate validation:

#### **Audio Files**
- **`data/synth/audio/lung_good.wav`**: High-quality lung audio (20s, 8kHz, SNR 24dB)
- **`data/synth/audio/heart_nearmiss.wav`**: Poor-quality heart audio (8s, 4kHz, SNR 12dB)

#### **Image Metadata**
- **`data/synth/images/retina_good.json`**: Good quality retina metadata
- **`data/synth/images/derm_nearmiss.json`**: Poor quality dermatology metadata

### **Usage**
```bash
# Generate synthetic test data
python3 scripts/make_synthetic_fixtures.py

# Test with readiness checkers
python3 tests/fot_image/tools/image_readiness_checker.py tests/fot_image/fixtures/retina_good.json
python3 tests/fot_audio/tools/audio_readiness_checker.py tests/fot_audio/fixtures/lung_good.json
```

---

## 📊 **Real Medical Datasets**

### **Available Datasets**

#### **🫀 Heart Sounds (PCG)**
- **Source**: PhysioNet/Computing in Cardiology Challenge 2016
- **Records**: 3,125 heart sound recordings
- **License**: Public training set
- **Usage**: Normal/abnormal classification, quality assessment
- **Access**: https://www.physionet.org/content/challenge-2016/1.0.0/

#### **🫁 Chest X-rays**
- **Source**: MIMIC-CXR Database
- **Records**: 377,110 chest X-ray images
- **License**: HIPAA de-identified, credentialed access required
- **Usage**: Chest pathology detection, image quality validation
- **Access**: https://physionet.org/content/mimic-cxr/

#### **🦠 Dermatoscopic Images**
- **Source**: HAM10000 Dataset
- **Records**: 10,015 skin lesion images
- **License**: Non-commercial use
- **Usage**: Skin lesion classification, dermoscopy analysis
- **Access**: https://complexity.cecs.ucf.edu/ham10000/

#### **👁️ Retina Fundus**
- **Source**: DRIVE Database
- **Records**: 40 fundus images with vessel ground truth
- **License**: Research/education use (no redistribution)
- **Usage**: Vessel segmentation, diabetic retinopathy screening
- **Access**: https://www.isi.uu.nl/research/databases/

#### **🫁 Respiratory Sounds**
- **Source**: ICBHI 2017 Respiratory Sound Database
- **Records**: 920 samples, 5.5 hours with crackle/wheeze labels
- **License**: Login required, widely cited
- **Usage**: Lung sound analysis, adventitious sound detection
- **Access**: https://ai4eu.dei.uc.pt/respiratory-sounds-dataset/

---

## 🔧 **Usage Instructions**

### **Image Readiness Testing**
```bash
# Test good quality image
python3 tests/fot_image/tools/image_readiness_checker.py tests/fot_image/fixtures/retina_good.json

# Test poor quality image
python3 tests/fot_image/tools/image_readiness_checker.py tests/fot_image/fixtures/derm_nearmiss.json

# Test custom image metadata
python3 tests/fot_image/tools/image_readiness_checker.py your_image_metadata.json
```

### **Audio Readiness Testing**
```bash
# Test good quality audio
python3 tests/fot_audio/tools/audio_readiness_checker.py tests/fot_audio/fixtures/lung_good.json

# Test poor quality audio
python3 tests/fot_audio/tools/audio_readiness_checker.py tests/fot_audio/fixtures/heart_nearmiss.json

# Test custom audio metadata
python3 tests/fot_audio/tools/audio_readiness_checker.py your_audio_metadata.json
```

### **Synthetic Data Generation**
```bash
# Generate all synthetic test data
python3 scripts/make_synthetic_fixtures.py

# View generated files
ls -la data/synth/
```

### **Real Dataset Download**
```bash
# View dataset information and download instructions
./scripts/get_test_data.sh

# Follow individual dataset instructions for download
```

---

## 📈 **Quality Standards**

### **Image Quality Thresholds**
- **Focus Score**: >0.6 (sharp, clear image)
- **Exposure Score**: 0.3-0.9 (proper lighting)
- **SNR**: >20 dB (good signal quality)
- **Resolution**: >512px shortest side (sufficient detail)

### **Audio Quality Thresholds**
- **SNR**: >20 dB (clear audio signal)
- **Noise Floor**: <-35 dBFS (low background noise)
- **Artifact Score**: <0.4 (minimal distortion)
- **Sample Rate**: >4000 Hz (adequate frequency range)
- **Duration**: >10 seconds (sufficient analysis time)

### **Technical Requirements**
- **PHI Compliance**: No patient identifiers
- **Calibration**: Device properly calibrated
- **Metadata**: Complete technical specifications
- **Format**: Standard medical file formats

---

## 🚀 **Integration with Clinical Pipeline**

### **Quality Gates**
The readiness checkers serve as quality gates in the clinical pipeline:

1. **Data Input**: Medical images/audio uploaded
2. **Quality Validation**: Readiness checker validates quality
3. **Gate Decision**: Pass/block based on quality metrics
4. **Quantum Processing**: Only high-quality data processed
5. **Clinical Output**: Reliable diagnostic recommendations

### **Automated Validation**
```python
# Example integration
from tests.fot_image.tools.image_readiness_checker import readiness_report

# Validate image before processing
image_metadata = load_image_metadata(image_file)
validation_result = readiness_report(image_metadata)

if validation_result["ready"]:
    # Process with quantum clinical engine
    quantum_result = process_with_quantum_engine(image_file)
else:
    # Request reacquisition or flag for review
    handle_poor_quality_data(validation_result["warnings"])
```

---

## 📚 **Documentation & Resources**

### **Technical Documentation**
- **Image Readiness Checker**: `tests/fot_image/tools/image_readiness_checker.py`
- **Audio Readiness Checker**: `tests/fot_audio/tools/audio_readiness_checker.py`
- **Synthetic Data Generator**: `scripts/make_synthetic_fixtures.py`
- **Test Datasets Guide**: `docs/TEST_DATASETS.md`

### **Test Fixtures**
- **Good Quality Examples**: Pass all quality gates
- **Near-Miss Examples**: Fail specific quality criteria
- **Synthetic Data**: Safe for demos and testing
- **Real Datasets**: Professional medical data

### **Quality Standards**
- **Medical Imaging**: Professional radiology standards
- **Medical Audio**: Clinical auscultation standards
- **PHI Compliance**: HIPAA privacy requirements
- **Technical Specifications**: Medical device standards

---

**🖼️🎵 Comprehensive Medical Data Quality Validation**

*Ensuring professional medical standards for images and audio in quantum clinical AI.*

---

**© 2024 Fortress AI - Medical Data Quality Validation Framework**
