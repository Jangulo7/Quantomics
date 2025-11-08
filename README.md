# Quantomics: Multi-Agent Protein Folding System

<p align="center">
  <img src="https://img.shields.io/badge/Hackathon-Qiskit%20Fall%20Fest%202025-blueviolet" alt="Hackathon">
  <img src="https://img.shields.io/badge/Event-Quantum%20Madrid%20×%20UPM-orange" alt="Event">
  <img src="https://img.shields.io/badge/Team-Quantomics-green" alt="Team">
  <img src="https://img.shields.io/badge/Status-Proof%20of%20Concept-yellow" alt="Status">
</p>

---

## 🎯 Project Overview

**Quantomics** presents a **multi-agent system architecture** for protein structure prediction that intelligently combines classical AI, quantum computing, and experimental validation. This project was developed for the **Qiskit Fall Fest Madrid 2025 - Quantum Madrid × UPM** hackathon.

### 🏗️ Project Status

- ✅ **Architecture Designed**: Complete multi-agent system architecture with 6 specialized agents
- ✅ **Individual Components Implemented**: Core functionality of each agent demonstrated
- ⏳ **Full Integration Pending**: Agent coordination and orchestration not yet implemented

> **Note**: This repository demonstrates a **proof of concept** with working implementations of individual agents. The complete multi-agent orchestration system represents future work.

---

## 🧬 Scientific Motivation

### The Protein Folding Challenge

Protein structure prediction remains one of biology's grand challenges. While recent AI models like AlphaFold2 and ESMFold have achieved remarkable success, they struggle with:

- **Low confidence regions** (intrinsically disordered regions)
- **Novel protein sequences** (orphan proteins lacking homologs)
- **Flexible loops** and dynamic regions
- **Variant-specific mutations** (e.g., SARS-CoV-2 variants)

### Our Innovation: Hybrid Classical-Quantum Approach

We propose an **intelligent routing system** that:

1. **Uses classical AI (ESMFold)** for fast, reliable predictions
2. **Routes to quantum computing (VQE)** when confidence is low (15-50%)
3. **Validates with experimental methods** (RDKit molecular analysis)

**Key Insight**: Not all proteins require quantum computing—use it only when classical methods show low confidence.

---

## 🏛️ System Architecture

### Multi-Agent System Design

```
User Input (Protein Sequence)
        ↓
[Coordination Agent] ← Orchestrates workflow
        ↓
[Verification Agent] ← Validates against UniProt/PDB
        ↓
[Classical AI Agent] ← ESMFold prediction
        ↓
    Confidence Analysis
        ├─ High (≥50%) → [Reporting Agent]
        ├─ Medium (15-50%) → [Quantum Agent] → [Reporting Agent]
        └─ Low (<15%) → Warning → Optional [Quantum Agent]
        ↓
[Experimental Agent] ← RDKit validation (if needed)
        ↓
    Final Report
```

### Agent Specifications

| Agent | Purpose | Status | Implementation |
|-------|---------|--------|----------------|
| **1. Coordination Agent** | Workflow orchestration & routing | 🟡 Designed | Architecture only |
| **2. Verification Agent** | Sequence validation (UniProt/PDB) | 🟡 Designed | Architecture only |
| **3. Classical AI Agent** | ESMFold structure prediction | ✅ **Implemented** | `esmfold.ipynb` |
| **4. Quantum Agent** | VQE quantum folding analysis | ✅ **Implemented** | `vqe_quiskit.py`, `quantum_protein_vqe.py` |
| **5. Reporting Agent** | Results aggregation & visualization | 🟡 Designed | Architecture only |
| **6. Experimental Agent** | RDKit molecular validation | ✅ **Implemented** | `Grafo.ipynb` |

### Decision Logic

The system uses **confidence-based routing** to optimize computational resources:

```python
if confidence >= 0.50:
    # HIGH confidence - Classical AI is sufficient
    return classical_result
    
elif 0.15 <= confidence < 0.50:
    # MEDIUM confidence - Quantum can help
    quantum_result = run_quantum_analysis()
    return combine(classical_result, quantum_result)
    
else:  # confidence < 0.15
    # LOW confidence - Warn user
    if user_accepts_risk:
        return run_quantum_analysis()  # May not help
    else:
        return classical_result
```

---

## 🔬 Implemented Components

### 1. Classical AI Agent (ESMFold)

**File**: `esmfold.ipynb`

**Implementation**:
- Uses Meta's ESMFold model for fast structure prediction
- Analyzes confidence scores (pLDDT) per residue
- Identifies high vs low confidence regions

**Test Cases**:
- ✅ **High Confidence**: SARS-CoV-2 RBD Core (223 residues)
  - Mean pLDDT: **74.5%**
  - Well-structured β-sheet domain
  - Range: 67-79%

- ⚠️ **Low Confidence**: SARS-CoV-2 Omicron BA.5 Flexible Loop (10 residues)
  - Mean pLDDT: **22.4%**
  - Intrinsically disordered region
  - Range: 16-30%

**Results Visualization**:

![ESMFold Confidence Comparison](esmfold_confidence_comparison.png)

**Key Findings**:
- ESMFold clearly distinguishes structured vs disordered regions
- Low pLDDT regions are ideal candidates for quantum analysis
- Δ pLDDT = 52.1 points between high and low confidence regions

---

### 2. Quantum Agent (VQE)

**Files**: 
- `quantum_protein_vqe.py` - Main implementation
- `vqe_quiskit.py` - Alternative implementation

**Method**: Variational Quantum Eigensolver (VQE) on 2D HP Lattice Model

**HP Lattice Model**:
```python
# Convert amino acids to Hydrophobic (H) or Polar (P)
Hydrophobic (H): A, V, I, L, M, F, Y, W
Polar (P):       R, N, D, C, Q, E, G, H, K, P, S, T

# Energy function
E = Σ ε_ij for H-H contacts
where ε_ij = -1.0 (favorable interaction)
```

**VQE Implementation**:
- **Quantum Circuit**: Hardware-efficient ansatz with RY rotations + CNOT gates
- **Optimizer**: COBYLA (200 iterations)
- **Backend**: Qiskit StatevectorEstimator (Qiskit 2.x compatible)
- **Runtime**: ~8 minutes for both regions

**Quantum Analysis Results**:

![VQE Quantum Analysis](vqe_protein_analysis.png)

**Results Summary**:

| Region | Ground State Energy | Stability Score | Interpretation |
|--------|---------------------|-----------------|----------------|
| **RBD Core** | -0.600 | 0.60 | LOW stability (unexpected) |
| **Flexible Loop** | -6.257 | 6.26 | HIGH stability (unexpected) |

**Biological Interpretation**:
⚠️ Unexpected results suggest:
- HP model is oversimplified for this application
- Need for more sophisticated Hamiltonians
- Importance of longer-range interactions beyond nearest neighbors
- Potential for improvement with real quantum hardware

**HP Patterns**:
```
RBD Core:        PHPPPPPH (8 residues analyzed)
Flexible Loop:   HPPHPHPHHP (10 residues)
```

**Key Insight**: The flexible loop has more hydrophobic residues (6 H vs 4 P), which in the simplified HP model creates more favorable interactions. This demonstrates the need for more realistic energy models that account for:
- Secondary structure propensity
- Long-range interactions
- Solvent accessibility
- Entropic contributions

---

### 3. Experimental Agent (RDKit)

**File**: `Grafo.ipynb`

**Implementation**:
- Molecular graph generation from SMILES codes
- Node and edge analysis (atoms, bonds)
- Distance matrix calculation
- Integration with VQE for energy optimization

**Capabilities**:
```python
class MolLoader:
    def load_molecules_QAOA()  # Load target molecules
    def load_smi()             # Parse SMILES format
    # Generates RDKit molecule objects for analysis
```

**Use Cases**:
1. **SMILES-based VQE Analysis**:
   - Convert protein to SMILES representation
   - Generate molecular graph
   - Extract atomic coordinates
   - Use distances in VQE Hamiltonian

2. **Graph Analysis**:
   - Node features (atom types, rings)
   - Edge features (bond types, distances)
   - Connectivity patterns

**Status**: Basic infrastructure implemented, integration with protein sequences pending

---

## 📊 Scientific Results

### ESMFold Confidence Analysis

**High Confidence Target** (SARS-CoV-2 RBD Core):
```
✅ Mean pLDDT: 74.5 ± 2.1%
✅ Range: 67-79%
✅ Confidence: HIGH
✅ Low confidence regions: 0 residues
✅ Interpretation: Well-folded, stable structure
```

**Low Confidence Target** (Omicron BA.5 Flexible Loop):
```
⚠️ Mean pLDDT: 22.4 ± 4.7%
⚠️ Range: 16-30%
⚠️ Confidence: LOW
⚠️ Low confidence regions: 10 residues (100%)
⚠️ Interpretation: Disordered, flexible region
```

### VQE Quantum Results

**Computational Details**:
- Qubits used: 8 (RBD Core fragment), 10 (Flexible Loop)
- Optimizer: COBYLA (200 iterations)
- Runtime: ~4-5 minutes per region
- Backend: Qiskit StatevectorEstimator

**Energy Landscape**:
```
RBD Core:        E = -0.600, Stability = 0.60 (LOW)
Flexible Loop:   E = -6.257, Stability = 6.26 (HIGH)
```

**Analysis**: 
The counterintuitive results (flexible loop showing higher quantum stability) highlight the limitations of simplified HP models and suggest areas for improvement:
1. Include backbone constraints
2. Add secondary structure penalties
3. Incorporate solvent effects
4. Use more sophisticated force fields

---

## 🚀 Getting Started

### Prerequisites

```bash
# Python 3.10+
python --version

# Required libraries
pip install qiskit qiskit-aer qiskit-algorithms
pip install torch fair-esm biotite
pip install rdkit matplotlib numpy pandas
```

### Running Individual Components

#### 1. ESMFold Analysis

```bash
# Open Jupyter notebook
jupyter notebook esmfold.ipynb

# Run all cells to:
# - Predict structures for test sequences
# - Calculate confidence scores
# - Visualize pLDDT distributions
```

#### 2. Quantum VQE Analysis

```bash
# Run VQE simulation (~8 minutes)
python quantum_protein_vqe.py

# Output:
# - Ground state energies
# - Stability scores
# - HP lattice patterns
# - Visualization: vqe_protein_analysis.png
```

#### 3. RDKit Molecular Analysis

```bash
# Open Jupyter notebook
jupyter notebook Grafo.ipynb

# Features:
# - Load SMILES molecules
# - Generate molecular graphs
# - Analyze atomic coordinates
```

---

## 📁 Repository Structure

```
quantomics/
├── README_HACKATHON.md                    # This file
├── architecture/
│   ├── protein_folding_architecture.mermaid   # Complete architecture diagram
│   ├── multiagent_architecture.py             # Full agent specifications
│   └── deployment_guide.md                    # Implementation roadmap
├── implementations/
│   ├── esmfold.ipynb                      # ✅ Classical AI Agent
│   ├── quantum_protein_vqe.py             # ✅ Quantum Agent (VQE)
│   ├── vqe_quiskit.py                     # ✅ Alternative VQE implementation
│   └── Grafo.ipynb                        # ✅ Experimental Agent (RDKit)
├── results/
│   ├── esmfold_confidence_comparison.png  # ESMFold results
│   └── vqe_protein_analysis.png           # VQE quantum results
└── docs/
    ├── ARCHITECTURE_SUMMARY.md            # Complete architecture documentation
    ├── implementation_examples.py         # Code examples for all agents
    └── sequence_diagram.md                # Message flow visualization
```

---

## 🔍 Key Findings & Insights

### ✅ What Works

1. **ESMFold Confidence Scoring**: 
   - Accurately identifies structured vs disordered regions
   - Clear separation: 74.5% vs 22.4% pLDDT
   - Reliable indicator for when to use quantum methods

2. **VQE Implementation**:
   - Successfully runs on Qiskit 2.x
   - Completes optimization in ~8 minutes
   - Demonstrates quantum energy calculations for proteins

3. **Modular Architecture**:
   - Each agent has well-defined responsibilities
   - Clean separation of concerns
   - Easy to extend or replace individual components

### ⚠️ Challenges & Limitations

1. **HP Model Oversimplification**:
   - Results don't match biological expectations
   - Need for more sophisticated energy models
   - Missing: backbone constraints, secondary structure, solvent effects

2. **Quantum-Classical Integration**:
   - Coordination between agents not implemented
   - Message passing infrastructure needed
   - State management required

3. **Scalability**:
   - VQE limited to ~10 qubits (10 residues)
   - Full proteins require fragmentation approach
   - Classical simulation is computationally expensive

### 💡 Future Improvements

#### Short-term (Next Hackathon)
- [ ] Implement agent coordination layer (Redis Pub/Sub)
- [ ] Add verification agent (UniProt/PDB queries)
- [ ] Improve HP model with secondary structure terms
- [ ] Create unified reporting system

#### Medium-term (6 months)
- [ ] Implement Miyazawa-Jernigan energy matrix
- [ ] Add RDKit SMILES-to-VQE pipeline
- [ ] Test on real IBM Quantum hardware
- [ ] Benchmark against AlphaFold2

#### Long-term (1 year)
- [ ] Replace HP model with full atomistic Hamiltonian
- [ ] Add molecular dynamics validation
- [ ] Scale to 20+ qubit proteins
- [ ] Publish scientific paper on hybrid approach

---

## 🎓 Educational Value

This project demonstrates several important concepts:

### Computer Science
- **Multi-agent architectures**: Coordination, communication, decision logic
- **Message-driven systems**: Pub/Sub patterns, event-driven design
- **Microservices**: Decoupled, scalable components

### Quantum Computing
- **VQE algorithm**: Variational quantum algorithms for optimization
- **Hamiltonian design**: Encoding physical problems for quantum computers
- **Hybrid quantum-classical**: Combining strengths of both paradigms

### Bioinformatics
- **Protein structure prediction**: Modern AI approaches (ESMFold)
- **Confidence metrics**: pLDDT scoring and interpretation
- **HP lattice model**: Simplified protein folding models

### Scientific Computing
- **Workflow orchestration**: Complex multi-step pipelines
- **Error handling**: Routing based on confidence/quality metrics
- **Validation**: Experimental verification of computational predictions

---

## 👥 Team Quantomics

**Hackathon**: Qiskit Fall Fest Madrid 2025 - Quantum Madrid × UPM

**Team Members**: Johanna Angulo, Enrique José Gómez García, Diego Barroso Martín, Adrián Madrid Gámez, Guillermo Carrasco Soto 

---

## 📚 References

### Classical AI
- **ESMFold**: Lin et al. (2023) "Evolutionary-scale prediction of atomic-level protein structure with a language model." Science.
- **AlphaFold2**: Jumper et al. (2021) "Highly accurate protein structure prediction with AlphaFold." Nature.

### Quantum Computing
- **VQE**: Peruzzo et al. (2014) "A variational eigenvalue solver on a photonic quantum processor." Nature Communications.
- **Protein Folding**: Perdomo-Ortiz et al. (2012) "Finding low-energy conformations of lattice protein models by quantum annealing." Scientific Reports.

### Protein Folding Models
- **HP Model**: Dill (1985) "Theory for the folding and stability of globular proteins." Biochemistry.
- **Miyazawa-Jernigan**: Miyazawa & Jernigan (1996) "Residue–residue potentials with a favorable contact pair term and an unfavorable high packing density term." J. Mol. Biol.

### Multi-Agent Systems
- **Agent Architectures**: Wooldridge (2009) "An Introduction to MultiAgent Systems." Wiley.
- **LangGraph**: LangChain AI (2024) "Building Stateful Multi-Actor Applications with LLMs."

---

## 📄 License

This project is licensed under the MIT License - see LICENSE file for details.

---

## 🙏 Acknowledgments

- **Quantum Madrid × UPM** for organizing Qiskit Fall Fest 2025
- **Meta AI** for the ESMFold model
- **IBM Quantum** for Qiskit framework
- **RDKit** community for molecular analysis tools
- **UniProt** and **PDB** for structural biology databases

---

## 🚀 Future Vision

This project represents the **first step** toward a production-ready hybrid classical-quantum system for protein structure prediction. We envision:

1. **Near-term**: Completing the multi-agent orchestration
2. **Mid-term**: Testing on real quantum hardware (IBM Quantum, IonQ)
3. **Long-term**: Deployment as a cloud service for the scientific community

**Join us in building the future of computational biology! 🧬⚛️**

---

<p align="center">
  <b>Developed with ❤️ by Team Quantomics</b><br>
  Qiskit Fall Fest Madrid 2025 - Quantum Madrid × UPM
</p>
