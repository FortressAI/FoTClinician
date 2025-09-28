# ☁️ FoTChemistry Streamlit Cloud Deployment Guide

## 🚀 Quick Deployment to Streamlit Cloud

### Step 1: Access Streamlit Cloud
1. Visit **https://share.streamlit.io**
2. Sign in with your GitHub account
3. Click **"New app"**

### Step 2: Configure App
- **Repository**: `FortressAI/FoTChemistry`
- **Branch**: `main`
- **Main file path**: `streamlit_app.py`
- **Advanced settings**: Set Python version to `3.9`

### Step 3: Requirements (Automatic)
The app will automatically use:
- **Main requirements**: `requirements-cloud.txt` (cloud-optimized)
- **System packages**: `packages.txt` (for RDKit dependencies)
- **Streamlit config**: `.streamlit/config.toml`

## 🔧 How Cloud Deployment Works

### Automatic Environment Detection
```python
# The app detects cloud deployment automatically:
is_cloud_deployment = (
    os.environ.get('STREAMLIT_SHARING') == '1' or
    os.environ.get('STREAMLIT_CLOUD') == '1' or
    'streamlit.io' in os.environ.get('HOSTNAME', '') or
    not os.path.exists('akg/client.py')  # No local AKG
)
```

### Data Loading Strategy
1. **Cloud**: Uses static data snapshot (`cloud_data_snapshot.json`)
2. **Local**: Uses live Neo4j data and file exports
3. **Fallback**: Embedded demo data (always works)

## 📊 Available Data

### Cloud Snapshot Features
- **43 unique molecules** from local discovery engine
- **High-quality discoveries** (score range: 0.577-0.786)
- **Complete molecular properties** (formula, weight, LogP, TPSA)
- **2D/3D visualization** (when RDKit available)

### Demo Molecules (Fallback)
- **Propanol** (CCCO) - Score: 0.786
- **Biphenyl** (c1ccc(-c2ccccc2)cc1) - Score: 0.754

## 🧬 Cloud App Features

### ✅ What Works in Cloud
- **Discovery dashboard** with all molecules
- **Molecular analysis** with property calculation
- **2D structures** (if RDKit installs successfully)
- **3D visualization** (if py3Dmol loads)
- **Property statistics** and scoring
- **Search and filtering**

### ⚠️ Cloud Limitations
- **No live discovery** (uses static snapshot)
- **No Neo4j connection** (not needed for display)
- **RDKit may not install** (graceful fallback to text)
- **Limited to snapshot data** (43 molecules)

## 🔄 Updating Cloud Data

### To Update the Cloud Snapshot:
1. **Run discovery locally** to generate new molecules
2. **Export new snapshot**: `python3 fix_molecular_generation.py`
3. **Commit and push**: Updates `cloud_data_snapshot.json`
4. **Streamlit Cloud** auto-deploys the new data

## 🎯 Expected Cloud URL
Your deployed app will be available at:
```
https://share.streamlit.io/fortressai/fotchemistry/main/streamlit_app.py
```

## 📞 Troubleshooting

### If RDKit Fails to Install
- App will show **text-based molecular data**
- **Fallback to SMILES strings** for structure display
- **All other features** remain functional

### If 3D Visualization Fails
- App will show **2D structures only**
- **Property analysis** still works
- **Molecular data** remains accessible

### If All Packages Fail
- App will use **embedded demo data**
- **Basic functionality** guaranteed
- **Contact support** if this happens

## ✅ Success Indicators

When deployment works correctly, you'll see:
- **"☁️ Cloud deployment detected"** message
- **"📦 Loaded 43 molecules from cloud_data_snapshot.json"**
- **Dashboard with molecular structures**
- **Working dropdown molecule selector**

---

**Ready to deploy? Visit https://share.streamlit.io and connect the FoTChemistry repository!** 🚀
