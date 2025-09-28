# 🧬 FoTChemistry User Guide

**Complete guide to running autonomous chemistry discovery locally**

## 🎯 Overview

FoTChemistry is a **perpetual discovery machine** that autonomously hunts, hypothesizes, tests, and collapses chemistry claims to truth. It connects to your existing Neo4j instance and provides real-time discovery visualization.

## ✅ Prerequisites

### Required
- **Neo4j running locally** (you already have this!)
  - Default: `bolt://localhost:7687` 
  - Username: `neo4j`
  - Password: `fotquantum`
- **Python 3.9+**
- **Basic dependencies**: `pip install streamlit plotly pandas neo4j pyyaml`

### Verified Working
✅ Neo4j connection tested and working  
✅ FoTChem schema safely isolated from protein/fluid apps  
✅ Discovery export pipeline functional  
✅ Streamlit dashboard operational  

---

## 🚀 Quick Start

### 1. **Test Your Connection**
```bash
python3 run_discovery_local.py --test-connection
```

Expected output:
```
🏥 AKG Health: Neo4j=True
✅ Neo4j connection successful!
📊 Current discoveries: X
```

### 2. **Start the Dashboard**
```bash
python3 -m streamlit run streamlit_app.py
```

Navigate to: **http://localhost:8501**

### 3. **Run Discovery (Export Only)**
```bash
python3 run_discovery_local.py --export-only
```

This exports existing data for the Streamlit dashboard without running new discovery.

---

## 🔬 Running Discovery Campaigns

### **Single Discovery Session**
```bash
# Run for 1 hour (default)
python3 run_discovery_local.py

# Run for 30 minutes
python3 run_discovery_local.py --budget 30m

# Run for 2 hours
python3 run_discovery_local.py --budget 2h
```

### **Continuous Discovery**
```bash
# Run continuously (restarts after each cycle)
python3 run_discovery_local.py --continuous --budget 1h
```

### **Specific Campaigns**
```bash
# Run only PFAS campaign
python3 run_discovery_local.py --campaigns configs/campaigns/pfas_leads.yaml

# Run multiple campaigns
python3 run_discovery_local.py --campaigns configs/campaigns/pfas_leads.yaml configs/campaigns/co2_catalysts.yaml
```

---

## 📊 Dashboard Features

### **Real-Time Metrics**
- 🏆 **Total Discoveries**: Claims collapsed to truth
- ✅ **Truth Collapsed**: Validated discoveries  
- ❌ **Refuted Claims**: Claims proven false
- ⏳ **Needs Evidence**: Claims requiring more data

### **Discovery Timeline**
- Chronological view of all discoveries
- Color-coded by verdict (truth/refuted/pending)
- Grouped by campaign

### **Connection Status**
- ✅ Neo4j connection indicator
- 🔄 Manual refresh controls
- 📁 Data export triggers

---

## 🗂️ Data Flow

### **Discovery Process**
1. **Scout** → Monitors data sources for new chemistry signals
2. **Architect** → Generates testable claims from signals  
3. **Alchemist** → Proposes candidate molecules/reactions
4. **Physicist** → Runs computational measurements
5. **Statistician** → Evaluates evidence vs. collapse rules
6. **Ethicist** → Screens for safety/ethical concerns
7. **Adjudicator** → Collapses claims to truth/refute/pending

### **Data Storage**
- **Neo4j**: Live discovery data with `FoTChem_` namespace
- **JSON Export**: `results/chemistry_discoveries.json` for Streamlit
- **Logs**: `logs/discovery_YYYYMMDD_HHMMSS.log`

---

## ⚙️ Configuration

### **Campaign Files**
Location: `configs/campaigns/*.yaml`

Example structure:
```yaml
name: "PFAS Sorbent Discovery"
objective: "minimize_residual_pfas_ngL"
batch_size: 16
virtue_weighting:
  Beneficence: 0.6
  Prudence: 0.3
  Honesty: 0.1
sources:
  - type: "kg-delta"
    query: "NEW_PFAS_STRUCTURES"
collapse_rules:
  success_if:
    residual_pfas_ngL: {"<=": 10}
  uncertainty_max: 0.2
ethics:
  blocklists: ["hazard_class_1"]
  allowlists: ["activated_carbons"]
```

### **Neo4j Configuration**
Edit `akg/client.py` if needed:
```python
'neo4j': {
    'uri': 'bolt://localhost:7687',
    'user': 'neo4j', 
    'password': 'fotquantum'  # Your password
}
```

---

## 🔧 Command Reference

### **Discovery Runner (`run_discovery_local.py`)**

| Command | Description |
|---------|-------------|
| `--test-connection` | Test Neo4j connection only |
| `--export-only` | Export existing data without discovery |
| `--continuous` | Run discovery continuously |
| `--budget 30m\|1h\|2h\|4h\|8h` | Set time budget |
| `--campaigns PATTERN` | Specify campaign files |

### **Streamlit Dashboard (`streamlit_app.py`)**

| Feature | Description |
|---------|-------------|
| Auto-refresh | Updates every 60 seconds |
| Manual refresh | Click "🔄 Refresh" button |
| Export trigger | Click "📁 Export Data" |
| Connection status | Real-time Neo4j health |

---

## 🛡️ Safety Features

### **Namespace Isolation**
- Uses `FoTChem_` prefixes for all chemistry nodes
- **SAFE** for existing protein/fluid apps
- No conflicts with your other FoT systems

### **Ethics Screening**
- Automatic screening of all generated candidates
- Configurable blocklists/allowlists
- Safety-first approach to discovery

### **Graceful Degradation**
- Works with Neo4j-only (GraphDB/Fuseki optional)
- Handles connection failures gracefully
- Continues from where it left off

---

## 📈 Monitoring Progress

### **Live Dashboard**
- **URL**: http://localhost:8501
- **Updates**: Real-time from Neo4j
- **Export**: Automatic JSON dumps for Git tracking

### **Log Files**
```bash
# View latest discovery log
tail -f logs/discovery_*.log

# Check for errors
grep ERROR logs/discovery_*.log
```

### **Neo4j Browser**
```bash
# View chemistry discoveries directly
MATCH (d:FoTChem_Discovery) RETURN d LIMIT 10

# Check claim status
MATCH (c:FoTChem_Claim) 
RETURN c.status, count(c) as count
```

---

## 🚨 Troubleshooting

### **Connection Issues**
```bash
# Test Neo4j connection
python3 run_discovery_local.py --test-connection

# Check if Neo4j is running
sudo netstat -tlnp | grep 7687
```

### **Discovery Not Starting**
```bash
# Check campaign files exist
ls configs/campaigns/*.yaml

# Verify agent modules
python3 -c "from agents.physicist.measure import MeasureQueue; print('OK')"
```

### **Streamlit Not Loading**
```bash
# Check if port is available
lsof -i :8501

# Start with specific port
streamlit run streamlit_app.py --server.port 8502
```

### **Data Not Updating**
```bash
# Force data export
python3 run_discovery_local.py --export-only

# Check JSON export
cat results/chemistry_discoveries.json | jq '.total_discoveries'
```

---

## 🎉 Success Indicators

### **Discovery Working**
- ✅ New entries in `logs/discovery_*.log`
- ✅ Increasing discovery count in dashboard
- ✅ New `FoTChem_` nodes in Neo4j
- ✅ Updated `results/chemistry_discoveries.json`

### **Dashboard Working**  
- ✅ Loads at http://localhost:8501
- ✅ Shows "✅ Neo4j Connected" in sidebar
- ✅ Displays discovery metrics > 0
- ✅ Refresh button updates data

### **System Health**
```bash
# Quick health check
python3 run_discovery_local.py --test-connection
streamlit --version
python3 -c "from akg.client import AKG; print('AKG OK')"
```

---

## 📞 Next Steps

1. **Start with test connection** to verify setup
2. **Run export-only** to populate dashboard
3. **Launch Streamlit** to view existing data
4. **Run short discovery** (30m) to test pipeline
5. **Monitor dashboard** for new discoveries
6. **Scale up** to continuous discovery

**You're now ready to run autonomous chemistry discovery! 🧬✨**
