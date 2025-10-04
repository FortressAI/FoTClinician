# 🔧 **FoTClinician Installation Guide**

## 🎯 **Installation Overview**

FoTClinician can be deployed via multiple methods depending on your infrastructure needs:

- **🌐 Streamlit Cloud** (Recommended for quick start)
- **🐳 Docker** (Containerized deployment)
- **🗄️ Local Installation** (Development/testing)
- **☁️ Enterprise Cloud** (AWS, Azure, GCP)

---

## 🌐 **Method 1: Streamlit Cloud Deployment (Recommended)**

### **⚡ Quick Start (5 minutes)**
1. **Create Streamlit Cloud Account**:
   - Go to [share.streamlit.io](https://share.streamlit.io)
   - Sign in with GitHub account
   - Authorize repository access

2. **Configure Deployment**:
   - Click "New app"
   - Repository: `FortressAI/FoTClinician`
   - Main file: `streamlit_app.py`
   - Python version: `3.9`

3. **Set Environment Variables**:
   ```toml
   # In Streamlit Cloud secrets
   [clinical]
   quantum_dimension = 512
   virtue_threshold = 0.25
   usmle_mode = validation
   ```

4. **Deploy**:
   - Click "Deploy!"
   - Wait 5-10 minutes
   - Access: `https://your-app.streamlit.app`

### **📊 Recommended Settings**
- **Memory**: 2GB
- **Storage**: 100MB
- **Timeout**: 60 seconds
- **Auto-restart**: Enabled

---

## 🐳 **Method 2: Docker Deployment**

### **📦 Using Docker Compose**

**Create `docker-compose.yml`:**

```yaml
version: '3.8'

services:
  fotclinician:
    image: fotclinician:latest
    build:
      context: .
      dockerfile: Dockerfile
    ports:
      - "8501:8501"
    environment:
      - VQBIT_DIMENSION=512
      - VIRTUE_THRESHOLD=0.25
      - USMLE_MODE=validation
    volumes:
      - ./data:/app/data
      - ./logs:/app/logs
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8501/_stcore/health"]
      interval: 30s
      timeout: 10s
      retries: 3

  redis:
    image: redis:7-alpine
    ports:
      - "6379:6379"
    volumes:
      - redis_data:/data
    restart: unless-stopped

volumes:
  redis_data:
```

**Deploy Application:**

```bash
# Clone repository
git clone https://github.com/FortressAI/FoTClinician.git
cd FoTClinician

# Build and start services
docker-compose up -d

# Check status
docker-compose ps
docker-compose logs -f fotclinician
```

### **🔧 Building Custom Image**

**Create `Dockerfile`:**

```dockerfile
FROM python:3.9-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Create non-root user
RUN useradd -m -u 1001 clinician
RUN chown -R clinician:clinician /app
USER clinician

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Expose port
EXPOSE 8501

# Start application
CMD ["streamlit", "run", "streamlit_app.py", "--server.host", "0.0.0.0", "--server.port", "8501"]
```

---

## 🗄️ **Method 3: Local Installation**

### **📋 Prerequisites**

**System Requirements:**
- Python 3.9+
- 2GB RAM minimum
- 500MB disk space
- Internet connection

**Operating Systems:**
- macOS 10.14+
- Ubuntu 18.04+
- Windows 10+

### **🐍 Python Installation**

**macOS (Using Homebrew):**
```bash
# Install Python 3.9
brew install python@3.9

# Verify installation
python3.9 --version

# Install pip
curl https://bootstrap.pypa.io/get-pip.py | python3.9 -
```

**Ubuntu/Debian:**
```bash
# Update package list
sudo apt update

# Install Python 3.9 and pip
sudo apt install python3.9 python3.9-pip python3.9-venv

# Verify installation
python3.9 --version
pip3 --version
```

**Windows:**
1. Download Python 3.9 from [python.org](https://www.python.org/downloads/)
2. Run installer with "Add to PATH" checked
3. Verify: `python --version` in Command Prompt

### **📦 Application Setup**

```bash
# Clone repository
git clone https://github.com/FortressAI/FoTClinician.git
cd FoTClinician

# Create virtual environment
python3.9 -m venv fotclinician_env

# Activate virtual environment

# macOS/Linux:
source fotclinician_env/bin/activate

# Windows:
fotclinician_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Create configuration
cp .streamlit/config.toml.template .streamlit/config.toml

# Start application
streamlit run streamlit_app.py
```

### **⚙️ Local Configuration**

**File**: `.streamlit/config.toml`

```toml
[global]
gatherUsageStats = false
showWarningOnDirectExecution = false

[server]
headless = false
port = 8501
enableCORS = true
enableXsrfProtection = true
maxUploadSize = 200

[browser]
gatherUsageStats = false
serverAddress = "localhost"

[theme]
primaryColor = "#1f77b4"
backgroundColor = "#ffffff"
secondaryBackgroundColor = "#f0f2f6"
textColor = "#262730"
```

---

## ☁️ **Method 4: Enterprise Cloud Deployment**

### **🌩️ AWS Deployment**

**EC2 Instance Setup:**
```bash
# Launch Ubuntu 20.04 instance (t3.medium recommended)
aws ec2 run-instances \
  --image-id ami-0c02fb55956c7d316 \
  --instance-type t3.medium \
  --key-name your-key-pair \
  --security-groups fotclinician-sg \
  --user-data file://install.sh
```

**Installation Script (`install.sh`):**
```bash
#!/bin/bash
yum update -y
yum install -y python39 git docker

# Clone repository
git clone https://github.com/FortressAI/FoTClinician.git
cd FoTClinician

# Install dependencies
pip3 install -r requirements.txt

# Create systemd service
cat > /etc/systemd/system/fotclinician.service << EOF
[Unit]
Description=FoTClinician Quantum Clinical Assistant
After=network.target

[Service]
Type=simple
User=ec2-user
WorkingDirectory=/home/ec2-user/FoTClinician
ExecStart=/usr/bin/python3.9 -m streamlit run streamlit_app.py --server.host 0.0.0.0 --server.port 8501
Restart=always

[Install]
WantedBy=multi-user.target
EOF

systemctl enable fotclinician
systemctl start fotclinician
```

### **🔵 Azure Deployment**

**Virtual Machine Setup:**
```bash
# Create resource group
az group create --name fotclinician-rg --location eastus

# Create VM
az vm create \
  --resource-group fotclinician-rg \
  --name fotclinician-vm \
  --image UbuntuLTS \
  --size Standard_D2s_v3 \
  --admin-username azureuser \
  --ssh-key-name your-key

# Install application
ssh azureuser@<vm-ip> "
  sudo apt update
  sudo apt install -y python3.9 python3.9-pip git
  git clone https://github.com/FortressAI/FoTClinician.git
  cd FoTClinician
  pip3 install -r requirements.txt
  nohup python3.9 -m streamlit run streamlit_app.py --server.host 0.0.0.0 --server.port 8501 &
"
```

### **🌐 Google Cloud Deployment**

**Compute Engine Setup:**
```bash
# Create VM instance
gcloud compute instances create fotclinician-vm \
  --zone=us-central1-a \
  --machine-type=e2-standard-2 \
  --image-family=ubuntu-2004-lts \
  --image-project=ubuntu-os-cloud

# Deploy application
gcloud compute ssh fotclinician-vm --command="
  sudo apt update && sudo apt install -y python3.9 python3.9-pip git
  git clone https://github.com/FortressAI/FoTClinician.git
  cd FoTClinician
  pip3 install -r requirements.txt
  nohup python3.9 -m streamlit run streamlit_app.py --server.host 0.0.0.0 --server.port 8501 &
"
```

---

## 🔧 **Configuration Options**

### **🔧 Environment Variables**

**Core Settings:**
```bash
# Quantum engine configuration
export VQBIT_DIMENSION=512
export QUANTUM_DECOHERENCE_RATE=0.1
export ENTANGLEMENT_STRATEGY=gradient

# Clinical AI settings
export USMLE_VALIDATION_MODE=strict
export SAFETY_PROTOCOL_LEVEL=maximum
export PASSING_THRESHOLD=0.85

# Performance settings
export MAX_CONCURRENT_USERS=100
export REQUEST_TIMEOUT=30
export CACHE_ENABLED=true
```

**Security Settings:**
```bash
# Authentication
export SECRET_KEY=your-secret-key-here
export JWT_SECRET=jwt-secret-key
export SESSION_TIMEOUT=3600

# HTTPS/Security
export HTTPS_ENABLED=true
export SSL_CERT_PATH=/path/to/cert.pem
export SSL_KEY_PATH=/path/to/key.pem
```

### **📊 Advanced Configuration**

**Quantum Engine Tuning:**
```yaml
# quantum_config.yaml
quantum_engine:
  vqbit_dimension: 512
  decoherence_rate: 0.1
  entanglement_threshold: 0.25
  measurement_fidelity: 0.99
  
virtue_supervisor:
  enabled: true
  hierarchy: ["non_maleficence", "honesty", "prudence", "justice"]
  threshold_override: false
  
clinical_ai:
  usmle_mode: validation
  specialty_weights:
    emergency_medicine: 1.0
    internal_medicine: 1.0
    pediatrics: 1.0
```

---

## 🔍 **Installation Verification**

### **🧪 Health Checks**

**Application Health:**
```bash
# Test local deployment
curl -f http://localhost:8501/_stcore/health

# Expected response: "ok"

# Test quantum engine
python3.9 -c "
from core.clinical.quantum_clinical_engine import QuantumClinicalEngine
engine = QuantumClinicalEngine(vqbit_dimension=256)
print('✅ Quantum engine working')
"

# Test data validator
python3.9 -c "
from core.clinical.data_readiness_checker import ClinicalDataContractValidator
checker = ClinicalDataContractValidator()
print('✅ Data validator working')
"
```

**Performance Benchmarks:**
```bash
# Start time test
time streamlit run streamlit_app.py --server.headless true --server.port 8510 &
sleep 10
curl http://localhost:8510/_stcore/health
kill %1

# Memory usage test
python3.9 -c "
import psutil
import sys
sys.path.append('.')
from core.clinical.quantum_clinical_engine import QuantumClinicalEngine
engine = QuantumClinicalEngine(vqbit_dimension=512)
print(f'Memory usage: {psutil.virtual_memory().percent}%')
"
```

### **🎯 Deployment Checklist**

- [ ] Application starts successfully
- [ ] Health endpoint responds with "ok"
- [ ] All imports working without errors
- [ ] Quantum engine initializes properly
- [ ] Data validator functions correctly
- [ ] Streamlit interface loads in browser
- [ ] Sample clinical case processes successfully
- [ ] USMLE validation tests pass
- [ ] Virtue supervision active
- [ ] Security settings enabled

---

## 🚨 **Troubleshooting**

### **❌ Common Issues**

**Import Errors:**
```bash
# If missing dependencies:
pip install --upgrade -r requirements.txt

# If Python version issues:
python3.9 --version  # Should be 3.9+
```

**Memory Issues:**
```bash
# Reduce quantum_dimension in config
echo 'VQBIT_DIMENSION=256' >> .streamlit/secrets.toml

# Increase swap space (Linux)
sudo fallocate -l 2G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

**Port Conflicts:**
```bash
# Find process using port 8501
lsof -i :8501

# Kill process and restart
kill -9 <PID>
streamlit run streamlit_app.py --server.port 8501
```

**Performance Issues:**
```bash
# Enable caching
echo 'GLOBAL_ENABLE_CACHING=true' >> .streamlit/secrets.toml

# Reduce log verbosity
echo 'STREAMLIT_LOGGER_LEVEL=WARNING' >> .streamlit/secrets.toml
```
