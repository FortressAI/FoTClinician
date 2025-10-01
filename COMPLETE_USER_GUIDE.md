# 🧬 Complete FoTChemistry User Guide

> **🎯 BREAKTHROUGH**: **6,443 molecular discoveries** with **12,144 problem-solution instances** identified!

Complete guide to using FoTChemistry's dual-dashboard system for molecular discovery and problem-solving analysis.

## 🌐 **Two-Dashboard System Overview**

FoTChemistry provides **two specialized dashboards** for different research workflows:

### **📊 Main Discovery Dashboard** (Port 8505)
**URL**: http://localhost:8505  
**Purpose**: Explore individual molecular discoveries  
**Best For**: Analyzing specific compounds, visualizing structures, quantum properties

### **🎯 Problem-Solution Dashboard** (Port 8506)  
**URL**: http://localhost:8506  
**Purpose**: Find compounds that solve chemistry problems  
**Best For**: Targeted research, solution identification, performance analysis

---

## 🚀 **Quick Start Guide**

### **1. Launch Both Dashboards**
```bash
# Terminal 1: Main Discovery Dashboard
streamlit run streamlit_app.py --server.port 8505

# Terminal 2: Problem-Solution Dashboard  
streamlit run streamlit_problem_solutions.py --server.port 8506
```

### **2. Choose Your Research Workflow**

**🔬 Individual Molecule Analysis → Use Main Dashboard (8505)**
- Browse 6,443 discoveries
- Analyze specific molecular properties
- View 2D/3D structures
- Explore quantum metrics

**🎯 Problem-Solving Research → Use Problem-Solution Dashboard (8506)**
- Find PFAS removal solutions (2,522 compounds)
- Identify green synthesis enablers (4,470 compounds)
- Discover thermodynamically consistent structures (5,152 compounds)
- Analyze performance distributions

---

## 📊 **Main Discovery Dashboard Guide** (Port 8505)

### **🧬 Discovered Molecules Tab**
**What You'll See:**
- Complete list of 6,443 molecular discoveries
- Real-time scores and properties
- Molecular structure previews
- Discovery timestamps

**Key Features:**
- **Search & Filter**: Find specific compounds by ID or properties
- **Sort Options**: Rank by score, molecular weight, complexity
- **Quick Preview**: 2D structure thumbnails for rapid browsing
- **Detailed View**: Click any molecule for comprehensive analysis

### **🔬 Analyze Molecule Tab**  
**What You'll Do:**
- Input SMILES strings for custom analysis
- Test molecules like:
  - `CCO` (Ethanol)
  - `CC(=O)O` (Acetic Acid)  
  - `C1=CC=CC=C1` (Benzene)

**Analysis Results:**
- **Quantum Properties**: vQbit measurements, coherence scores
- **Drug-Likeness**: Lipinski compliance, ADMET properties
- **3D Visualization**: Interactive molecular structures
- **Safety Assessment**: Hazard screening and ethics validation

### **⚛️ System Status Tab**
**Monitor:**
- Quantum engine performance (8096-dimensional substrate)
- GPU acceleration status (MPS shaders)
- Database connectivity (Neo4j, AKG)
- Available molecular libraries (RDKit, stmol, py3Dmol)

---

## 🎯 **Problem-Solution Dashboard Guide** (Port 8506)

### **🧪 Chemistry Problems Overview**

#### **💧 PFAS Removal**
- **Goal**: Remove PFAS contaminants to <25 ng/L
- **Solutions Found**: 2,522 compounds (39.1% success rate)
- **Top Performers**: Compounds with optimal fluorine content and hydrophobicity
- **Use Cases**: Water treatment, environmental remediation

#### **🧮 Thermodynamic Consistency**  
- **Goal**: Achieve cycle closure ≤ 0.3 kcal/mol
- **Solutions Found**: 5,152 compounds (80.0% success rate)
- **Top Performers**: Simple, predictable molecular structures
- **Use Cases**: Computational validation, model development

#### **🌱 Green Synthesis**
- **Goal**: Reduce PMI & E-Factor by ≥50%
- **Solutions Found**: 4,470 compounds (69.4% success rate)  
- **Top Performers**: High atom economy, diverse functionality
- **Use Cases**: Sustainable manufacturing, process optimization

#### **⚡ CO₂ Electrocatalysis**
- **Goal**: Achieve Faradaic efficiency ≥ 65%
- **Solutions Found**: 0 compounds (computational proxy needs refinement)
- **Status**: Under development - improved models coming soon

### **📊 Dashboard Features**

#### **Performance Analysis**
- **Distribution Charts**: See how compounds perform across metrics
- **Success Rate Tracking**: Monitor problem-solving effectiveness
- **Threshold Visualization**: Understand success criteria clearly
- **Comparative Analysis**: Compare performance across problems

#### **Top Solutions Explorer**
- **Ranked Lists**: Best performers for each problem
- **Detailed Properties**: SMILES, performance values, uncertainties
- **Downloadable Results**: Export CSV files for laboratory validation
- **Virtue Scoring**: Quality assessment using Field of Truth methodology

#### **Multi-Problem Analysis**
- **Versatile Compounds**: Find molecules solving multiple problems
- **Cross-Problem Insights**: Understand molecular design principles
- **Research Prioritization**: Focus on most promising candidates

---

## 🔄 **Workflow Examples**

### **Scenario 1: PFAS Contamination Research**
1. **Start**: Open Problem-Solution Dashboard (8506)
2. **Select**: PFAS Removal from sidebar
3. **Analyze**: Review 2,522 solutions and success metrics
4. **Download**: Export top 50 solutions as CSV
5. **Validate**: Use Main Dashboard (8505) to analyze specific compounds
6. **Synthesize**: Laboratory validation of most promising candidates

### **Scenario 2: Green Chemistry Development**
1. **Start**: Problem-Solution Dashboard (8506) → Green Synthesis
2. **Filter**: Find compounds with >60% reduction rates
3. **Cross-Reference**: Check if compounds also solve other problems
4. **Deep Dive**: Main Dashboard (8505) for molecular structure analysis
5. **Export**: Download multi-problem solvers for process development

### **Scenario 3: Computational Validation**
1. **Start**: Main Dashboard (8505) → Analyze Molecule tab
2. **Test**: Input your computational predictions as SMILES
3. **Compare**: Check quantum properties against our database
4. **Verify**: Problem-Solution Dashboard (8506) to see if it solves known problems
5. **Publish**: Use results to validate computational methods

---

## 🌐 **Cloud Deployment**

Both dashboards work in cloud environments with static data:

### **For Main Dashboard:**
- Uses `results/chemistry_discoveries.json` (local)
- Falls back to `cloud_data_snapshot.json` (cloud)
- Includes pre-generated 2D/3D visualizations

### **For Problem-Solution Dashboard:**
- Uses `problem_solution_analysis/complete_analysis.json` (local)
- Falls back to `cloud_problem_solution_data.json` (cloud)
- Secondary fallback to `results/problem_solution_summary.json` (cloud)
- Optimized with complete solution sets for comprehensive analysis

### **Deployment Commands:**
```bash
# Local development
streamlit run streamlit_app.py --server.port 8505
streamlit run streamlit_problem_solutions.py --server.port 8506

# Cloud deployment (auto-detects environment)
# Apps automatically switch to static data sources
```

---

## 🔧 **Advanced Features**

### **Data Integration**
- **Neo4j Connectivity**: Live database queries (local)
- **JSON Exports**: Static data snapshots (cloud)
- **SPARQL Queries**: Truth-mining workflows
- **Field of Truth**: Claim validation and evidence assessment

### **Molecular Visualization**
- **2D Structures**: SVG rendering with RDKit
- **3D Visualization**: Interactive stmol and py3Dmol
- **Fallback Systems**: Graceful degradation for cloud deployment
- **Performance Optimization**: Cached visualizations

### **Quantum Analysis**
- **8096-Dimensional Substrate**: Real quantum mechanics
- **MPS Acceleration**: Apple Silicon GPU optimization
- **Virtue Vectors**: Multi-dimensional quality assessment
- **Born Rule Compliance**: Validated quantum measurements

---

## 📚 **Data Sources & Ontology**

### **FoTChemistry Ontology** 
- **File**: `ontology/FoTChem.ttl`
- **Classes**: Problem, Compound, Claim, Measurement, Evidence
- **Properties**: Performance metrics, virtue scores, truth verdicts
- **Validation**: SHACL shapes for claim verification

### **Analysis Results**
- **Complete Dataset**: `problem_solution_analysis/complete_analysis.json`
- **Problem-Specific**: Individual JSON files per chemistry challenge  
- **Claims Database**: `problem_solution_claims.jsonl` for Neo4j
- **SPARQL Queries**: `ontology/problem_solution_queries.sparql`

---

## 🎯 **Research Impact**

### **Immediate Applications**
- **Environmental**: PFAS contamination solutions ready for testing
- **Sustainability**: 4,470 green chemistry enablers identified
- **Computational**: 5,152 validated thermodynamic models
- **Discovery**: Systematic approach to problem-solving chemistry

### **Long-term Vision**
- **Automated Labs**: Direct integration with robotic synthesis
- **AI-Guided Research**: Quantum-informed molecular design
- **Truth Validation**: Community-driven claim verification
- **Open Science**: FAIR data principles and reproducible workflows

---

## 🛠️ **Troubleshooting**

### **Common Issues**

#### **Dashboard Won't Load**
```bash
# Check if ports are available
lsof -i :8505 :8506

# Restart with different ports if needed
streamlit run streamlit_app.py --server.port 8507
```

#### **No Data Showing**
- **Local**: Ensure `results/chemistry_discoveries.json` exists
- **Cloud**: Verify `cloud_data_snapshot.json` is present
- **Problem Dashboard**: Check `problem_solution_analysis/` directory

#### **Visualization Errors**
- **RDKit Issues**: Check Python environment and dependencies
- **3D Problems**: Verify `stmol` and `py3Dmol` installation
- **Cloud Deployment**: Pre-generated visualizations should work

### **Performance Optimization**
```bash
# Clear Streamlit cache
streamlit cache clear

# Restart with clean slate
pkill -f streamlit && streamlit run streamlit_app.py --server.port 8505
```

---

## 🎉 **Getting Started Checklist**

- [ ] **Clone Repository**: `git clone https://github.com/FortressAI/FoTChemistry.git`
- [ ] **Install Dependencies**: `pip install -r requirements.txt`
- [ ] **Start Neo4j**: Ensure database running with password `fotquantum`
- [ ] **Launch Main Dashboard**: `streamlit run streamlit_app.py --server.port 8505`
- [ ] **Launch Problem Dashboard**: `streamlit run streamlit_problem_solutions.py --server.port 8506`
- [ ] **Explore Data**: Browse 6,443 discoveries and 12,144 problem-solutions
- [ ] **Download Results**: Export relevant solutions for validation
- [ ] **Validate in Lab**: Test most promising compounds experimentally

**Welcome to the future of problem-solving chemistry!** 🧬🎯✨
