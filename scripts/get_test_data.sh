#!/usr/bin/env bash
# FoT Clinical Test Data Downloader
# Downloads small, legally redistributable medical datasets for testing

set -euo pipefail

echo "🔬 FoT Clinical Test Data Downloader"
echo "====================================="

# Create data directory
mkdir -p data/tests/{images,audio,synthetic}

echo ""
echo "📊 Available Test Datasets:"
echo ""

# Heart sounds (PCG) - PhysioNet Challenge 2016
echo "🫀 Heart Sounds (PCG):"
echo "   Source: PhysioNet/Computing in Cardiology Challenge 2016"
echo "   Records: 3,125 heart sound recordings"
echo "   License: Public training set"
echo "   Download: https://www.physionet.org/content/challenge-2016/1.0.0/"
echo "   Usage: Normal/abnormal classification, quality assessment"
echo ""

# Chest X-rays - MIMIC-CXR
echo "🫁 Chest X-rays (MIMIC-CXR):"
echo "   Source: MIMIC-CXR Database"
echo "   Records: Large HIPAA de-identified chest X-ray dataset"
echo "   License: Requires credentialed PhysioNet access"
echo "   Download: https://physionet.org/content/mimic-cxr/"
echo "   Usage: Chest pathology detection, image quality validation"
echo ""

# Dermatoscopic images - HAM10000
echo "🦠 Dermatoscopic Images (HAM10000):"
echo "   Source: Complex Adaptive Systems Laboratory"
echo "   Records: 10,015 skin lesion images"
echo "   License: Non-commercial use"
echo "   Download: https://complexity.cecs.ucf.edu/ham10000/"
echo "   Usage: Skin lesion classification, dermoscopy analysis"
echo ""

# Retina fundus - DRIVE
echo "👁️ Retina Fundus (DRIVE):"
echo "   Source: Image Sciences Institute, Utrecht"
echo "   Records: 40 fundus images with vessel ground truth"
echo "   License: Research/education use (no redistribution)"
echo "   Download: https://www.isi.uu.nl/research/databases/"
echo "   Usage: Vessel segmentation, diabetic retinopathy screening"
echo ""

# Respiratory sounds - ICBHI 2017
echo "🫁 Respiratory Sounds (ICBHI 2017):"
echo "   Source: ICBHI Respiratory Sound Database"
echo "   Records: 920 samples, 5.5 hours with crackle/wheeze labels"
echo "   License: Login required, widely cited"
echo "   Download: https://ai4eu.dei.uc.pt/respiratory-sounds-dataset/"
echo "   Usage: Lung sound analysis, adventitious sound detection"
echo ""

echo "📋 Download Instructions:"
echo "========================="
echo ""
echo "1. For PhysioNet datasets:"
echo "   - Create account at https://physionet.org/"
echo "   - Complete CITI training if required"
echo "   - Download using wget or PhysioNet tools"
echo ""
echo "2. For HAM10000:"
echo "   - Use Dataverse or Kaggle CLI"
echo "   - Respect non-commercial license terms"
echo ""
echo "3. For DRIVE:"
echo "   - Manual download from ISI Utrecht"
echo "   - Use only for research/education"
echo ""
echo "4. For ICBHI:"
echo "   - Create account and download audio + annotations"
echo "   - Follow citation requirements"
echo ""

echo "⚠️  Important Notes:"
echo "==================="
echo ""
echo "- These datasets are for TESTING and VALIDATION only"
echo "- Do NOT commit actual medical data to repository"
echo "- Always respect dataset licenses and terms of use"
echo "- Use synthetic data for demos and CI/CD pipelines"
echo "- Real datasets should be downloaded locally by users"
echo ""

echo "🔧 Synthetic Data Generation:"
echo "============================="
echo ""
echo "Run the synthetic data generator for demo purposes:"
echo "python scripts/make_synthetic_fixtures.py"
echo ""

echo "✅ Test data downloader ready!"
echo "   Use individual dataset links above to download legally."
