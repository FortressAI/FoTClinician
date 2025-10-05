# 📊 **FoT Clinical Test Datasets**

**Comprehensive Guide to Medical Test Data for FoTClinician Validation**

---

## 🎯 **Overview**

This document provides detailed information about medical datasets suitable for testing and validating the FoTClinician quantum clinical pipeline. All datasets listed are either publicly available, HIPAA de-identified, or have clear licensing terms for research use.

---

## 🫀 **Heart Sounds (Phonocardiogram)**

### **PhysioNet/Computing in Cardiology Challenge 2016**

**📋 Dataset Details:**
- **Records**: 3,125 heart sound recordings
- **Duration**: Variable (typically 5-30 seconds)
- **Quality**: Mixed environments, real-world conditions
- **Labels**: Normal/abnormal classification
- **Format**: WAV files with metadata

**🔗 Access:**
- **URL**: https://www.physionet.org/content/challenge-2016/1.0.0/
- **License**: Public training set, research use
- **Requirements**: PhysioNet account (free)
- **Download**: Direct via PhysioNet tools

**🎯 Clinical Applications:**
- Heart murmur detection
- Cardiac abnormality screening
- Audio quality assessment
- Signal processing validation

**📊 Quality Metrics Available:**
- Signal-to-noise ratio
- Artifact detection
- Duration analysis
- Environmental noise assessment

---

## 🫁 **Chest X-rays**

### **MIMIC-CXR Database**

**📋 Dataset Details:**
- **Records**: 377,110 chest X-ray images
- **Format**: DICOM and JPG variants available
- **Labels**: Structured findings and impressions
- **De-identification**: HIPAA compliant
- **Quality**: Hospital-grade imaging

**🔗 Access:**
- **URL**: https://physionet.org/content/mimic-cxr/
- **License**: Requires credentialed PhysioNet access
- **Requirements**: CITI training completion
- **Download**: Via PhysioNet tools after approval

**🎯 Clinical Applications:**
- Chest pathology detection
- Image quality validation
- Radiologist agreement studies
- AI model benchmarking

**📊 Quality Metrics Available:**
- Image resolution analysis
- Exposure quality assessment
- Artifact detection
- Anatomical completeness

---

## 🦠 **Dermatoscopic Images**

### **HAM10000 Dataset**

**📋 Dataset Details:**
- **Records**: 10,015 skin lesion images
- **Categories**: 7 different lesion types
- **Resolution**: Variable (typically 600x450)
- **Format**: JPG images with metadata
- **Annotations**: Expert dermatologist labels

**🔗 Access:**
- **URL**: https://complexity.cecs.ucf.edu/ham10000/
- **License**: Non-commercial use
- **Requirements**: Registration required
- **Download**: Dataverse or Kaggle mirror

**🎯 Clinical Applications:**
- Skin lesion classification
- Melanoma detection
- Dermoscopy analysis
- Teledermatology validation

**📊 Quality Metrics Available:**
- Focus quality assessment
- Lighting analysis
- Color accuracy
- Lesion visibility

---

## 👁️ **Retina Fundus Images**

### **DRIVE Database**

**📋 Dataset Details:**
- **Records**: 40 fundus images
- **Resolution**: 565x584 pixels
- **Format**: TIFF images
- **Annotations**: Vessel segmentation ground truth
- **Quality**: Research-grade imaging

**🔗 Access:**
- **URL**: https://www.isi.uu.nl/research/databases/
- **License**: Research/education use (no redistribution)
- **Requirements**: Manual download
- **Download**: Direct from ISI Utrecht

**🎯 Clinical Applications:**
- Vessel segmentation
- Diabetic retinopathy screening
- Image quality assessment
- Ophthalmology AI validation

**📊 Quality Metrics Available:**
- Vessel visibility
- Image sharpness
- Color balance
- Field of view coverage

---

## 🫁 **Respiratory Sounds**

### **ICBHI 2017 Respiratory Sound Database**

**📋 Dataset Details:**
- **Records**: 920 audio samples
- **Duration**: 5.5 hours total
- **Labels**: Crackle/wheeze annotations
- **Quality**: Expert-labeled recordings
- **Format**: WAV files with annotations

**🔗 Access:**
- **URL**: https://ai4eu.dei.uc.pt/respiratory-sounds-dataset/
- **License**: Login required, widely cited
- **Requirements**: Account creation
- **Download**: Direct after registration

**🎯 Clinical Applications:**
- Adventitious sound detection
- Lung disease screening
- Audio quality assessment
- Respiratory monitoring

**📊 Quality Metrics Available:**
- Signal-to-noise ratio
- Artifact detection
- Duration analysis
- Environmental noise

---

## 🔧 **Synthetic Test Data**

### **FoT Synthetic Generator**

**📋 Generated Content:**
- **Audio**: Synthetic heart and lung sounds
- **Images**: Placeholder medical images
- **Metadata**: Realistic clinical annotations
- **Quality**: Controlled test scenarios

**🎯 Usage:**
- **Testing**: Pipeline validation
- **Demo**: Safe demonstration data
- **CI/CD**: Automated testing
- **Development**: Local development

**📊 Generated Files:**
```
data/synth/
├── audio/
│   ├── lung_good.wav
│   └── heart_nearmiss.wav
├── images/
│   ├── retina_good.json
│   └── derm_nearmiss.json
└── metadata/
    └── synthetic_data_manifest.json
```

---

## 🚀 **Quick Start Guide**

### **1. Download Test Data**

```bash
# Run the test data downloader
./scripts/get_test_data.sh

# Generate synthetic data
python scripts/make_synthetic_fixtures.py
```

### **2. Validate Data Quality**

```bash
# Test image readiness
python tests/fot_image/tools/image_readiness_checker.py tests/fot_image/fixtures/retina_good.json

# Test audio readiness  
python tests/fot_audio/tools/audio_readiness_checker.py tests/fot_audio/fixtures/lung_good.json
```

### **3. Run CI Tests**

```bash
# Execute smoke tests
python -m pytest tests/ -v

# Run GitHub Actions locally
act -j fot_tests
```

---

## ⚖️ **Legal & Ethical Considerations**

### **🔒 Data Privacy**
- **HIPAA Compliance**: All datasets are de-identified
- **No PHI**: No personally identifiable information
- **Research Use**: Appropriate for research and testing
- **Commercial Restrictions**: Some datasets have usage limitations

### **📋 License Compliance**
- **PhysioNet**: Requires account and training
- **HAM10000**: Non-commercial use only
- **DRIVE**: Research/education, no redistribution
- **ICBHI**: Login required, citation needed

### **🎯 Best Practices**
- **Local Storage**: Download datasets locally
- **No Commits**: Don't commit actual medical data
- **Synthetic First**: Use synthetic data for demos
- **Respect Terms**: Follow all dataset licenses

---

## 📈 **Performance Benchmarks**

### **Expected Readiness Scores**

**✅ Good Quality Examples:**
- **Retina Fundus**: Focus >0.8, SNR >25dB
- **Lung Audio**: SNR >20dB, Duration >10s
- **Heart Sounds**: Calibration passed, Low artifacts

**⚠️ Near-Miss Examples:**
- **Derm Images**: Focus <0.6, Overexposed
- **Heart Audio**: SNR <15dB, Short duration
- **Poor Quality**: High noise, Artifacts present

### **Validation Targets**
- **Image Readiness**: >90% pass rate for good examples
- **Audio Readiness**: >85% pass rate for quality recordings
- **PHI Detection**: 100% blocking of flagged content
- **Quality Gates**: Appropriate warnings for substandard data

---

## 🔄 **Continuous Integration**

### **Automated Testing**
- **Smoke Tests**: Basic functionality validation
- **Quality Gates**: Readiness checker validation
- **Regression Tests**: Prevent quality degradation
- **Performance Tests**: Response time validation

### **CI Pipeline**
```yaml
# .github/workflows/fot_tests.yml
name: FoT Clinical Smoke Tests
on: [push, pull_request]
jobs:
  smoke:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Run readiness checks
        run: |
          python tests/fot_image/tools/image_readiness_checker.py tests/fot_image/fixtures/retina_good.json
          python tests/fot_audio/tools/audio_readiness_checker.py tests/fot_audio/fixtures/lung_good.json
```

---

## 📞 **Support & Resources**

### **🔗 Dataset Links**
- [PhysioNet](https://physionet.org/) - Medical data repository
- [Complexity Lab](https://complexity.cecs.ucf.edu/) - HAM10000 dataset
- [ISI Utrecht](https://www.isi.uu.nl/) - DRIVE database
- [AI4EU](https://ai4eu.dei.uc.pt/) - ICBHI respiratory sounds

### **📚 Documentation**
- [PhysioNet Documentation](https://physionet.org/about/)
- [Dataset Citation Guidelines](https://physionet.org/citation/)
- [FoT Clinical Testing Guide](TESTING_GUIDE.md)

### **🆘 Troubleshooting**
- **Download Issues**: Check dataset access requirements
- **Quality Failures**: Review readiness checker output
- **License Questions**: Contact dataset maintainers
- **Technical Support**: GitHub Issues repository

---

**🎯 Ready to validate your FoTClinician pipeline with real medical data!**

*Use these datasets responsibly and in accordance with all licensing terms and ethical guidelines.*
