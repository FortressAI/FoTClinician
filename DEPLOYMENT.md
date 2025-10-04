# 🚀 FoTClinician Deployment Guide

## 📦 Professional Deployment Instructions for Quantum Medical AI

---

## 🎯 **System Requirements**

### Minimum Requirements

- **🖥️ Operating System**: Ubuntu 20.04+, Windows 10+, macOS 10.15+
- **🐍 Python**: 3.8 or higher
- **💾 RAM**: 8GB minimum, 16GB recommended
- **💿 Storage**: 5GB free space
- **🌐 Network**: Internet connection for models/dependencies

### Recommended Requirements  

- **🖥️ OS**: Ubuntu 22.04 LTS
- **🐍 Python**: 3.10+
- **💾 RAM**: 32GB
- **💿 Storage**: 50GB SSD
- **🌐 Network**: High-bandwidth connection

---

## 🛠️ **Installation Methods**

### Method 1: Quick Start (Development)

```bash
# Clone repository
git clone https://github.com/FortressAI/FoTClinician.git
cd FoTClinician

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Launch application
streamlit run streamlit_app.py
```

**🌐 Access**: Navigate to `http://localhost:8501`

### Method 2: Container Deployment (Production)

```bash
# Pull Docker image (if available)
docker pull fortressai/fotclinician:latest

# Run container
docker run -d \
  --name fotclinician \
  -p 8501:8501 \
  -e VIRTUE_SUPERVISOR_ENABLED=true \
  -e USMLE_VALIDATION_MODE=strict \
  fortressai/fotclinician:latest
```

### Method 3: Cloud Platform Deployment

#### 🚀 **AWS EC2 Deployment**

```bash
# Launch EC2 instance
aws ec2 run-instances \
  --image-id ami-0c02fb55956c7d316 \
  --instance-type t3.xlarge \
  --key-name your-key \
  --security-groups fotclinician-sg

# Connect and deploy
ssh -i your-key.pem ubuntu@your-ec2-ip
sudo apt update && sudo apt install -y python3-pip docker.io
git clone https://github.com/FortressAI/FoTClinician.git
cd FoTClinician
pip3 install -r requirements.txt
streamlit run streamlit_app.py --server.host 0.0.0.0
```

#### 🌐 **Google Cloud Platform**

```bash
# Create Compute Engine instance
gcloud compute instances create fotclinician-vm \
  --machine-type=e2-standard-4 \
  --image-family=ubuntu-2004-lts \
  --image-project=ubuntu-os-cloud

# Deploy application
gcloud compute ssh fotclinician-vm --command="
  sudo apt update && sudo apt install -y python3-pip git
  git clone https://github.com/FortressAI/FoTClinician.git
  cd FoTClinician
  pip3 install -r requirements.txt
  nohup streamlit run streamlit_app.py --server.host 0.0.0.0 --server.port 8501 &
"
```

#### 🔵 **Microsoft Azure**

```bash
# Create VM
az vm create \
  --resource-group fotclinician-rg \
  --name fotclinician-vm \
  --image UbuntuLTS \
  --size Standard_D4s_v3 \
  --admin-username azureuser \
  --ssh-key-name your-ssh-key

# Deploy via SSH
ssh azureuser@your-vm-ip "
  sudo apt update && sudo apt install -y python3-pip git
  git clone https://github.com/FortressAI/FoTClinician.git
  cd FoTClinician && pip3 install -r requirements.txt
  nohup streamlit run streamlit_app.py --server.host 0.0.0.0 &
"
```

---

## ⚙️ **Configuration Management**

### Environment Variables

```bash
# Core quantum engine settings
export VQBIT_DIMENSION=512
export QUANTUM_DECOHERENCE_RATE=0.1
export VIRTUE_SUPERVISOR_ENABLED=true

# Medical validation settings
export USMLE_VALIDATION_MODE=strict
export SAFETY_PROTOCOL_LEVEL=maximum
export PASSING_THRESHOLD=0.85

# Security settings
export SECRET_KEY=your-secret-key
export HTTPS_ENABLED=true
export SSL_CERT_PATH=/path/to/cert.pem
export SSL_KEY_PATH=/path/to/key.pem

# Database settings (if applicable)
export DATABASE_URL=postgresql://user:pass@localhost/fotclinician
export REDIS_URL=redis://localhost:6379
```

### Configuration Files

**File**: `/etc/fotclinician/config.yaml`

```yaml
# Quantum Clinical Engine Configuration
quantum_engine:
  vqbit_dimension: 512
  decoherence_rate: 0.1
  entanglement_strategy: "gradient"
  measurement_fidelity: 0.99

# Virtue Supervisor Settings  
virtue_supervisor:
  enabled: true
  hierarchy: ["non_maleficence", "honesty", "prudence", "justice"]
  conflict_resolution: "weighted_average"
  audit_logging: true

# Medical Validation
validation:
  usmle_mode: "strict"
  passing_threshold: 0.85
  safety_level: "maximum"
  comprehensive_testing: true

# Security Configuration
security:
  api_authentication: "token-based"
  request_rate_limiting: "100/minute"
  audit_trail: true
  encryption:
    algorithm: "AES-256"
    key_rotation_days: 30

# Logging Configuration  
logging:
  level: "INFO"
  format: "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
  handlers:
    - type: "file"
      path: "/var/log/fotclinician/app.log"
      max_bytes: 10485760
      backup_count: 5
    - type: "stream"
      stream: "stderr"
```

---

## 🔧 **Advanced Deployment**

### Multi-Instance Deployment

#### Load Balancer Configuration

**Nginx Configuration**: `/etc/nginx/sites-available/fotclinician`

```nginx
upstream fotclinician_backend {
    server 127.0.0.1:8501 weight=3;
    server 127.0.0.1:8502 weight=3;
    server 127.0.0.1:8503 weight=2;
}

server {
    listen 80;
    server_name fotclinician.example.com;
    
    location / {
        proxy_pass http://fotclinician_backend;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_read_timeout 300s;
        proxy_connect_timeout 75s;
    }
    
    # Health check endpoint
    location /health {
        access_log off;
        return 200 "healthy\n";
        add_header Content-Type text/plain;
    }
}
```

#### Docker Compose Configuration

**File**: `docker-compose.yml`

```yaml
version: '3.8'

services:
  fotclinician:
    build: .
    ports:
      - "8501:8501"
    environment:
      - VIRTUE_SUPERVISOR_ENABLED=true
      - USMLE_VALIDATION_MODE=strict
      - DATABASE_URL=postgresql://postgres:password@db:5432/fotclinician
    volumes:
      - ./logs:/app/logs
      - ./data:/app/data
    depends_on:
      - db
      - redis

  db:
    image: postgres:13
    environment:
      POSTGRES_DB: fotclinician
      POSTGRES_USER: postgres
      POSTGRES_PASSWORD: password
    volumes:
      - postgres_data:/var/lib/postgresql/data

  redis:
    image: redis:6-alpine
    volumes:
      - redis_data:/data

  nginx:
    image: nginx:alpine
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./nginx.conf:/etc/nginx/nginx.conf
      - ./ssl:/etc/nginx/ssl
    depends_on:
      - fotclinician

volumes:
  postgres_data:
  redis_data:
```

### Kubernetes Deployment

**File**: `k8s/fotclinician-deployment.yaml`

```yaml
apiVersion: apps/v1
kind: Deployment
metadata:
  name: fotclinician
  labels:
    app: fotclinician
spec:
  replicas: 3
  selector:
    matchLabels:
      app: fotclinician
  template:
    metadata:
      labels:
        app: fotclinician
    spec:
      containers:
      - name: fotclinician
        image: fortressai/fotclinician:latest
        ports:
        - containerPort: 8501
        env:
        - name: VIRTUE_SUPERVISOR_ENABLED
          value: "true"
        - name: USMLE_VALIDATION_MODE  
          value: "strict"
        resources:
          requests:
            memory: "2Gi"
            cpu: "1000m"
          limits:
            memory: "4Gi"
            cpu: "2000m"
        livenessProbe:
          httpGet:
            path: /health
            port: 8501
          initialDelaySeconds: 30
          periodSeconds: 10
        readinessProbe:
          httpGet:
            path: /health  
            port: 8501
          initialDelaySeconds: 5
          periodSeconds: 5

---
apiVersion: v1
kind: Service
metadata:
  name: fotclinician-service
spec:
  selector:
    app: fotclinician
  ports:
    - protocol: TCP
      port: 80
      targetPort: 8501
  type: LoadBalancer
```

---

## 📊 **Monitoring & Maintenance**

### Health Monitoring

**Health Check Script**: `scripts/health_check.py`

```python
#!/usr/bin/env python3
import requests
import json
import sys

def health_check():
    """Comprehensive health check for FoTClinician"""
    
    try:
        # Check Streamlit app
        response = requests.get('http://localhost:8501/_stcore/health', timeout=10)
        if response.status_code != 200:
            print("❌ Streamlit health check failed")
            return 1
        
        # Check quantum engine
        test_case = {
            'case_id': 'HEALTH_CHECK',
            'chief_complaint': 'test',
            'age': 30
        }
        
        # Import and test quantum engine
        from core.clinical.quantum_clinical_engine import QuantumClinicalEngine
        engine = QuantumClinicalEngine()
        quantum_case = engine.encode_clinical_case(test_case)
        
        if not quantum_case.differential_qbits:
            print("❌ Quantum engine health check failed")
            return 1
            
        print("✅ All health checks passed")
        return 0
        
    except Exception as e:
        print(f"❌ Health check failed: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(health_check())
```

### Automated Backup

**Backup Script**: `scripts/backup.py`

```python
#!/usr/bin/env python3
import os
import shutil
import datetime
import zipfile

def backup_fotclinician():
    """Create automated backup of FoTClinician"""
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_file = f"fotclinician_backup_{timestamp}.zip"
    
    with zipfile.ZipFile(backup_file, 'w') as backup:
        # Backup configuration files
        backup.write('config.yaml', 'config/config.yaml')
        
        # Backup logs
        if os.path.exists('logs/'):
            for root, dirs, files in os.walk('logs/'):
                for file in files:
                    backup.write(os.path.join(root, file))
        
        # Backup data
        if os.path.exists('data/'):
            for root, dirs, files in os.walk('data/'):
                for file in files:
                    backup.write(os.path.join(root, file))
    
    print(f"✅ Backup created: {backup_file}")
    return backup_file

if __name__ == "__main__":
    backup_fotclinician()
```

### Performance Monitoring

**Monitoring Script**: `scripts/monitor.py`

```python
#!/usr/bin/env python3
import psutil
import time
import json

def monitor_performance():
    """Monitor FoTClinician performance metrics"""
    
    metrics = {
        'timestamp': time.time(),
        'cpu_percent': psutil.cpu_percent(),
        'memory_percent': psutil.virtual_memory().percent,
        'disk_percent': psutil.disk_usage('/').percent,
        'processes': len(psutil.pids())
    }
    
    # Log metrics
    with open('performance_monitor.log', 'a') as f:
        f.write(json.dumps(metrics) + '\n')
    
    # Alert on high resource usage
    if metrics['cpu_percent'] > 80:
        print(f"⚠️ High CPU usage: {metrics['cpu_percent']:.1f}%")
    
    if metrics['memory_percent'] > 85:
        print(f"⚠️ High memory usage: {metrics['memory_percent']:.1f}%")
    
    return metrics

if __name__ == "__main__":
    while True:
        monitor_performance()
        time.sleep(60)  # Monitor every minute
```

---

## 🔒 **Security Configuration**

### SSL/TLS Setup

```bash
# Generate SSL certificate (self-signed for development)
openssl req -x509 -newkey rsa:4096 -keyout key.pem -out cert.pem -days 365 -nodes

# Launch with HTTPS
streamlit run streamlit_app.py \
  --server.sslCertFile=cert.pem \
  --server.sslKeyFile=key.pem \
  --server.port 8501
```

### Authentication Setup

**File**: `auth/middleware.py`

```python
from functools import wraps
import jwt
from flask import request, jsonify

def token_required(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        token = request.headers.get('Authorization')
        
        if not token:
            return jsonify({'message': 'Token is missing'}), 401
        
        try:
            data = jwt.decode(token, app.config['SECRET_KEY'], algorithms=['HS256'])
            current_user = data['user_id']
        except:
            return jsonify({'message': 'Token is invalid'}), 401
        
        return f(current_user, *args, **kwargs)
    
    return decorated
```

### Rate Limiting

```python
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

limiter = Limiter(
    app,
    key_func=get_remote_address,
    default_limits=["100 per minute"]
)

@app.route('/analyze')
@limiter.limit("10 per minute")  
def analyze_endpoint():
    # Quantum clinical analysis endpoint
    pass
```

---

## 📚 **Additional Resources**

### Support Contacts

- 🐛 **Bug Reports**: [GitHub Issues](https://github.com/FortressAI/FoTClinician/issues)
- 💬 **Community Support**: [GitHub Discussions](https://github.com/FortressAI/FoTClinician/discussions)
- 📧 **Enterprise Support**: fortress.ai@enterprise.support

### Documentation Links

- 📖 **[User Guide](docs/USER_GUIDE.md)**: Complete user manual
- 🔬 **[API Documentation](docs/API_DOCUMENTATION.md)**: Technical API reference
- 🛡️ **[Security Guide](docs/SECURITY_GUIDE.md)**: Security best practices
- ⚛️ **[Quantum Reference](docs/QUANTUM_REFERENCE.md)**: Scientific methodology

### Training Resources

- 🎓 **[USMLE Validation Guide](docs/USMLE_CERTIFICATION.md)**: Medical exam validation
- 🧪 **[Testing Framework](tests/)**: Comprehensive test suite
- 📊 **[Performance Optimization](docs/PERFORMANCE_GUIDE.md)**: Performance tuning

---

## ⚖️ **Legal & Compliance**

### Healthcare Compliance

- ✅ **HIPAA**: Configure appropriate data protection
- ✅ **FDA Guidance**: Follow AI/ML device regulations  
- ✅ **Medical Standards**: Maintain USMLE validation
- ✅ **Ethical AI**: Enforce virtue-based constraints

### License & Usage

**License**: MIT License - See [LICENSE](LICENSE) file for details
**Usage**: Educational and research use approved
**Commercial Use**: Contact Fortress AI for commercial licensing

---

*🚀 FoTClinician Professional Deployment Guide v1.0*

**© 2024 Fortress AI - Quantum Clinical Intelligence Platform**
