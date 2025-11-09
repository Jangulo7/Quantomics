"""
Quantum Protein Folding Analysis with VQE
==========================================
Analyzes SARS-CoV-2 protein regions using Variational Quantum Eigensolver
on a 2D HP lattice model. Runs in ~8 minutes on classical simulator.

Requirements:
pip install qiskit qiskit-aer qiskit-algorithms matplotlib numpy scipy
"""

"""
Quantum Protein Folding Analysis with VQE
==========================================
Analyzes SARS-CoV-2 protein regions using Variational Quantum Eigensolver
on a 2D HP lattice model. Runs in ~8 minutes on classical simulator.

FIXED for Qiskit 2.x compatibility
"""

import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit_aer import AerSimulator
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA, SLSQP

# FIXED: Import correct Estimator for Qiskit 2.x
from qiskit.primitives import StatevectorEstimator as Estimator

import time
import os

#

# ============================================================================
# PROTEIN DATA
# ============================================================================

HIGH_CONFIDENCE_TARGET = {
    'name': 'SARS-CoV-2 RBD Core (Structured)',
    'sequence': 'RVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF',
    'length': 223,
    'region': 'Full RBD domain',
}

LOW_CONFIDENCE_TARGET = {
    'name': 'SARS-CoV-2 Omicron BA.5 Flexible Loop',
    'sequence': 'YQPYRVVVLS',
    'length': 10,
    'region': 'Extended loop in RBD (residues 493-502)',
}

# ============================================================================
# HP MODEL: Convert amino acids to Hydrophobic (1) or Polar (0)
# ============================================================================

def sequence_to_hp(sequence):
    """
    Convert amino acid sequence to HP (Hydrophobic-Polar) model.
    
    Hydrophobic (H=1): A, V, I, L, M, F, Y, W
    Polar (P=0): R, N, D, C, Q, E, G, H, K, P, S, T
    """
    hydrophobic = set('AVILMFYW')
    hp_sequence = []
    
    for aa in sequence:
        if aa in hydrophobic:
            hp_sequence.append(1)  # H
        else:
            hp_sequence.append(0)  # P
    
    return np.array(hp_sequence)

def print_hp_sequence(name, sequence, hp_seq):
    """Visualize HP mapping"""
    print(f"\n{'='*80}")
    print(f"🧬 {name}")
    print(f"{'='*80}")
    print(f"Length: {len(sequence)} residues")
    print(f"\nHP Pattern:")
    
    # Print first 50 for visualization
    display_length = min(50, len(sequence))
    print(f"Sequence: {sequence[:display_length]}...")
    print(f"HP Code:  {''.join(['H' if h else 'P' for h in hp_seq[:display_length]])}...")
    
    h_count = np.sum(hp_seq)
    p_count = len(hp_seq) - h_count
    print(f"\n📊 Composition:")
    print(f"   Hydrophobic (H): {h_count} ({100*h_count/len(hp_seq):.1f}%)")
    print(f"   Polar (P):       {p_count} ({100*p_count/len(hp_seq):.1f}%)")
    
    return hp_seq

# ============================================================================
# 2D LATTICE HAMILTONIAN
# ============================================================================

def create_lattice_hamiltonian(hp_sequence, contact_energy=-1.0):
    """
    Create a simplified Hamiltonian for HP lattice model.
    
    Energy function:
    E = Σ ε_ij for all non-adjacent contacts
    where ε_ij = -1 if both H-H contact, 0 otherwise
    
    For simplicity, we'll create a toy model with n qubits representing
    n positions, and measure interactions between H residues.
    """
    n = len(hp_sequence)
    
    # For a tractable VQE problem, we'll encode interactions as Pauli strings
    # Each qubit represents whether a hydrophobic residue is "active" in forming contacts
    
    pauli_list = []
    
    # H-H contact interactions (simplified: nearest neighbors + diagonal)
    for i in range(n - 1):
        if hp_sequence[i] == 1 and hp_sequence[i+1] == 1:
            # H-H contact: favorable interaction
            # Z_i Z_{i+1} measures correlation
            pauli_str = ['I'] * n
            pauli_str[i] = 'Z'
            pauli_str[i+1] = 'Z'
            pauli_list.append((''.join(pauli_str), contact_energy))
    
    # Add some longer-range interactions for realism
    for i in range(n - 2):
        if hp_sequence[i] == 1 and hp_sequence[i+2] == 1:
            pauli_str = ['I'] * n
            pauli_str[i] = 'Z'
            pauli_str[i+2] = 'Z'
            pauli_list.append((''.join(pauli_str), contact_energy * 0.5))
    
    # Entropic penalty (disfavors too many active contacts - represents entropy loss)
    for i in range(n):
        if hp_sequence[i] == 1:
            pauli_str = ['I'] * n
            pauli_str[i] = 'Z'
            pauli_list.append((''.join(pauli_str), 0.3))
    
    if len(pauli_list) == 0:
        # No interactions (all polar) - add identity
        pauli_list.append(('I' * n, 0.0))
    
    hamiltonian = SparsePauliOp.from_list(pauli_list)
    
    return hamiltonian

# ============================================================================
# MJ MODEL: Consider the interaction energies of the Ma
# ============================================================================

#First, load the Miyazawa-Jernigan matrix as a tuple of ((Interaction energies)array, (aminoacid)list)

dir_path = os.path.dirname(os.path.realpath(__file__))

def load_mj_matrix(path):
    file_path = os.path.realpath(os.path.join(os.path.dirname(__file__), "mj_matrix.txt"))
    matrix = np.loadtxt(fname=file_path, dtype=str)
    energy_matrix = np.zeros((np.shape(matrix)[0], np.shape(matrix)[1]))
    for row in range(1, np.shape(matrix)[0]):
        for col in range(row - 1, np.shape(matrix)[1]):
            energy_matrix[row, col] = float(matrix[row, col])
    energy_matrix = energy_matrix[1:,]
    symbols = list(matrix[0, :])
    return energy_matrix, symbols

#After this, develop the Hamiltonian

def create_lattice_hamiltonian_mj(aa_sequence, energy_matrix, symbols, contact_energy=-1.0):

    n = len(aa_sequence)
    
    # For a tractable VQE problem, we'll encode interactions as Pauli strings
    # Each qubit represents whether a hydrophobic residue is "active" in forming contacts
    
    pauli_list = []
    
    # H-H contact interactions (simplified: nearest neighbors + diagonal)
    for i in range(n - 1):
        
        for j in range(n - 2):
            if i != j and i < j:
                pauli_str = ['I'] * n
                first_index = symbols.index(aa_sequence[i])
                second_index = symbols.index(aa_sequence[j+1])
                interaction_energy = energy_matrix[first_index, second_index]
                pauli_str[i] = 'Z'
                pauli_str[j] = 'Z'
                pauli_list.append((''.join(pauli_str), interaction_energy * (0.5 + 0.5 * (j-i) / n)))


    # Entropic penalty (disfavors too many active contacts - represents entropy loss)
    

    hamiltonian = SparsePauliOp.from_list(pauli_list)
    
    return hamiltonian

# ============================================================================
# VARIATIONAL ANSATZ
# ============================================================================

def create_ansatz(num_qubits, layers=2):
    """
    Create a hardware-efficient ansatz for VQE.
    Uses RY rotations and CNOT entanglement.
    """
    qc = QuantumCircuit(num_qubits)
    params = []
    
    # Initial layer of rotations
    for i in range(num_qubits):
        param = Parameter(f'θ_init_{i}')
        qc.ry(param, i)
        params.append(param)
    
    # Entangling layers
    for layer in range(layers):
        # Entanglement
        for i in range(num_qubits - 1):
            qc.cx(i, i + 1)
        
        # Rotations
        for i in range(num_qubits):
            param = Parameter(f'θ_L{layer}_{i}')
            qc.ry(param, i)
            params.append(param)
    
    return qc, params

# ============================================================================
# VQE OPTIMIZATION
# ============================================================================

def run_vqe_analysis(name, hp_sequence, max_qubits=8):
    """
    Run VQE to find ground state energy of HP lattice model.
    """
    print(f"\n{'='*80}")
    print(f"⚛️  QUANTUM VQE ANALYSIS: {name}")
    print(f"{'='*80}")
    
    # Truncate if too long
    if len(hp_sequence) > max_qubits:
        print(f"⚠️  Sequence too long ({len(hp_sequence)} residues)")
        print(f"   Using first {max_qubits} residues for quantum analysis")
        hp_sequence = hp_sequence[:max_qubits]
    
    num_qubits = len(hp_sequence)
    print(f"\n🔬 VQE Setup:")
    print(f"   Qubits: {num_qubits}")
    print(f"   HP Sequence: {''.join(['H' if h else 'P' for h in hp_sequence])}")
    
    # Create Hamiltonian
    print(f"\n📐 Building Hamiltonian...")
    hamiltonian = create_lattice_hamiltonian(hp_sequence)
    print(f"   Pauli terms: {len(hamiltonian)}")
    print(f"   Sample terms: {list(hamiltonian.to_list())[:3]}")
    
    # Create ansatz
    print(f"\n🎛️  Creating variational circuit...")
    ansatz, params = create_ansatz(num_qubits, layers=2)
    num_params = len(params)
    print(f"   Parameters: {num_params}")
    print(f"   Circuit depth: {ansatz.depth()}")
    
    # Setup VQE - FIXED for Qiskit 2.x
    print(f"\n⚡ Running VQE optimization...")
    estimator = Estimator()  # StatevectorEstimator for Qiskit 2.x
    optimizer = COBYLA(maxiter=200)
    
    # Initial point (small random values)
    np.random.seed(42)
    initial_point = np.random.uniform(-0.1, 0.1, num_params)
    
    vqe = VQE(estimator, ansatz, optimizer, initial_point=initial_point)
    
    start_time = time.time()
    result = vqe.compute_minimum_eigenvalue(hamiltonian)
    end_time = time.time()
    
    # Extract results
    energy = result.optimal_value
    optimal_params = result.optimal_point
    
    print(f"\n✅ VQE Complete!")
    print(f"   Time: {end_time - start_time:.2f} seconds")
    print(f"   Ground state energy: {energy:.4f}")
    print(f"   Optimizer iterations: {result.cost_function_evals}")
    
    # Analyze stability
    stability_score = -energy  # Lower energy = more stable = higher score
    
    print(f"\n📊 Structural Analysis:")
    print(f"   Stability score: {stability_score:.4f}")
    
    if stability_score > 2.0:
        stability = "HIGH (well-folded, stable)"
    elif stability_score > 1.0:
        stability = "MEDIUM (partially stable)"
    else:
        stability = "LOW (disordered, flexible)"
    
    print(f"   Interpretation: {stability}")
    
    return {
        'name': name,
        'hp_sequence': hp_sequence,
        'num_qubits': num_qubits,
        'energy': energy,
        'stability_score': stability_score,
        'stability': stability,
        'optimal_params': optimal_params,
        'time': end_time - start_time,
        'iterations': result.cost_function_evals
    }

def run_vqe_analysis_mj(name, aa_sequence, energy_matrix, symbols, max_qubits=8):
    """
    Run VQE to find ground state energy of HP lattice model.
    """
    print(f"\n{'='*80}")
    print(f"⚛️  QUANTUM VQE ANALYSIS: {name}")
    print(f"{'='*80}")
    
    # Truncate if too long
    if len(aa_sequence) > max_qubits:
        print(f"⚠️  Sequence too long ({len(aa_sequence)} residues)")
        print(f"   Using first {max_qubits} residues for quantum analysis")
        aa_sequence = aa_sequence[:max_qubits]
    
    num_qubits = len(aa_sequence)
    print(f"\n🔬 VQE Setup:")
    print(f"   Qubits: {num_qubits}")
    #print(f"   HP Sequence: {''.join(['H' if h else 'P' for h in aa_sequence])}")
    
    # Create Hamiltonian
    print(f"\n📐 Building Hamiltonian...")
    hamiltonian = create_lattice_hamiltonian_mj(aa_sequence, energy_matrix, symbols)
    print(f"   Pauli terms: {len(hamiltonian)}")
    print(f"   Sample terms: {list(hamiltonian.to_list())[:3]}")
    
    # Create ansatz
    print(f"\n🎛️  Creating variational circuit...")
    ansatz, params = create_ansatz(num_qubits, layers=2)
    num_params = len(params)
    print(f"   Parameters: {num_params}")
    print(f"   Circuit depth: {ansatz.depth()}")
    
    # Setup VQE - FIXED for Qiskit 2.x
    print(f"\n⚡ Running VQE optimization...")
    estimator = Estimator()  # StatevectorEstimator for Qiskit 2.x
    optimizer = COBYLA(maxiter=200)
    
    # Initial point (small random values)
    np.random.seed(42)
    initial_point = np.random.uniform(-0.1, 0.1, num_params)
    
    vqe = VQE(estimator, ansatz, optimizer, initial_point=initial_point)
    
    start_time = time.time()
    result = vqe.compute_minimum_eigenvalue(hamiltonian)
    end_time = time.time()
    
    # Extract results
    energy = result.optimal_value
    optimal_params = result.optimal_point
    
    print(f"\n✅ VQE Complete!")
    print(f"   Time: {end_time - start_time:.2f} seconds")
    print(f"   Ground state energy: {energy:.4f}")
    print(f"   Optimizer iterations: {result.cost_function_evals}")
    
    # Analyze stability
    stability_score = -energy  # Lower energy = more stable = higher score
    
    print(f"\n📊 Structural Analysis:")
    print(f"   Stability score: {stability_score:.4f}")
    
    if stability_score > 2.0:
        stability = "HIGH (well-folded, stable)"
    elif stability_score > 1.0:
        stability = "MEDIUM (partially stable)"
    else:
        stability = "LOW (disordered, flexible)"
    
    print(f"   Interpretation: {stability}")
    
    return {
        'name': name,
        'hp_sequence': aa_sequence,
        'num_qubits': num_qubits,
        'energy': energy,
        'stability_score': stability_score,
        'stability': stability,
        'optimal_params': optimal_params,
        'time': end_time - start_time,
        'iterations': result.cost_function_evals
    }


# ============================================================================
# VISUALIZATION
# ============================================================================



def plot_comparison(results_high, results_low):
    """
    Compare quantum VQE results for high vs low confidence regions.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('🧬 Quantum VQE (Miyazawa-Jernigan): SARS-CoV-2 Protein Regions', 
                 fontsize=16, fontweight='bold')
    
    # Stability comparison
    ax1 = axes[0, 0]
    names = [results_high['name'].split('(')[0].strip(), 
             results_low['name'].split('(')[0].strip()]
    stabilities = [results_high['stability_score'], results_low['stability_score']]
    colors = ['green' if s > 1.5 else 'orange' if s > 0.8 else 'red' for s in stabilities]
    
    bars = ax1.bar(names, stabilities, color=colors, alpha=0.7, edgecolor='black')
    ax1.set_ylabel('Stability Score (Higher = More Stable)', fontweight='bold')
    ax1.set_title('Ground State Stability Comparison')
    #ax1.axhline(y=1.5, color='green', linestyle='--', alpha=0.3, label='High stability threshold')
    #ax1.axhline(y=0.8, color='orange', linestyle='--', alpha=0.3, label='Medium stability threshold')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # Add value labels
    for bar, val in zip(bars, stabilities):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.2f}',
                ha='center', va='bottom', fontweight='bold')
    
    # Energy comparison
    ax2 = axes[0, 1]
    energies = [results_high['energy'], results_low['energy']]
    bars = ax2.bar(names, energies, color=['steelblue', 'coral'], alpha=0.7, edgecolor='black')
    ax2.set_ylabel('Ground State Energy (Lower = More Stable)', fontweight='bold')
    ax2.set_title('VQE Ground State Energy')
    ax2.grid(axis='y', alpha=0.3)
    
    for bar, val in zip(bars, energies):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.3f}',
                ha='center', va='bottom' if val < 0 else 'top', fontweight='bold')
    
    # FIXED: Amino acid sequence patterns (convert to numeric for visualization)
    ax3 = axes[1, 0]
    aa_high = results_high['hp_sequence']  # This is amino acid string like 'RVQPTESI'
    
    # Convert amino acid sequence to numeric HP values for visualization
    hp_high_numeric = sequence_to_hp(aa_high)
    
    ax3.imshow([hp_high_numeric], cmap='RdYlGn', aspect='auto', interpolation='nearest')
    ax3.set_yticks([0])
    ax3.set_yticklabels([results_high['name'].split('(')[0].strip()])
    ax3.set_title('HP Pattern: High Confidence Region')
    ax3.set_xticks(range(len(hp_high_numeric)))
    ax3.set_xticklabels(['H' if h else 'P' for h in hp_high_numeric], fontsize=8)
    ax3.set_xlabel('Residue Position')
    
    ax4 = axes[1, 1]
    aa_low = results_low['hp_sequence']  # This is amino acid string like 'YQPYRVVVLS'
    
    # Convert amino acid sequence to numeric HP values for visualization
    hp_low_numeric = sequence_to_hp(aa_low)
    
    ax4.imshow([hp_low_numeric], cmap='RdYlGn', aspect='auto', interpolation='nearest')
    ax4.set_yticks([0])
    ax4.set_yticklabels([results_low['name'].split('(')[0].strip()])
    ax4.set_xlabel('Residue Position')
    ax4.set_title('HP Pattern: Low Confidence Region')
    ax4.set_xticks(range(len(hp_low_numeric)))
    ax4.set_xticklabels(['H' if h else 'P' for h in hp_low_numeric], fontsize=10)
    
    plt.tight_layout()
    plt.savefig('vqe_protein_analysis.png', dpi=300, bbox_inches='tight')
    print("\n💾 Visualization saved as 'vqe_protein_analysis.png'")
    plt.show()



# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    print("\n" + "="*80)
    print("⚛️  QUANTUM PROTEIN FOLDING WITH VQE")
    print("="*80)
    print("\n🔬 Method: Variational Quantum Eigensolver (VQE)")
    print("📐 Model: HP Lattice (Hydrophobic-Polar) AND MJ interactions Lattice")
    print("💻 Backend: Qiskit StatevectorEstimator (Qiskit 2.x)")
    print("⏱️  Target runtime: ~8 minutes")
    
    # Convert sequences to HP model
    hp_high = print_hp_sequence(HIGH_CONFIDENCE_TARGET['name'],HIGH_CONFIDENCE_TARGET['sequence'],sequence_to_hp(HIGH_CONFIDENCE_TARGET['sequence']))
    
    hp_low = print_hp_sequence(LOW_CONFIDENCE_TARGET['name'], LOW_CONFIDENCE_TARGET['sequence'], sequence_to_hp(LOW_CONFIDENCE_TARGET['sequence']))
    
    print("\n" + "="*80)
    print("🚀 Starting VQE Quantum Simulations...")
    print("="*80)
    
    # Load MJ model
    energy_matrix, symbols = load_mj_matrix(dir_path)

    print("Energy matrix loaded!")

    #mj_high = create_lattice_hamiltonian_mj(HIGH_CONFIDENCE_TARGET['sequence'], energy_matrix, symbols, contact_energy=-1.0)
    #mj_low = create_lattice_hamiltonian_mj(LOW_CONFIDENCE_TARGET['sequence'], energy_matrix, symbols, contact_energy=-1.0)

    #print("MJ-High = " + str(mj_high))
    #print("MJ-Low = " + str(mj_low))

    # Run VQE for high confidence region (fragment)
    results_high = run_vqe_analysis(HIGH_CONFIDENCE_TARGET['name'], hp_high, max_qubits=8)  # Use 8-residue fragment
    print("Results high done!")
    results_high_mj = run_vqe_analysis_mj(HIGH_CONFIDENCE_TARGET['name'], HIGH_CONFIDENCE_TARGET['sequence'], energy_matrix, symbols, max_qubits=8)
    print("Results high - mj done!")
    # Run VQE for low confidence region
    results_low = run_vqe_analysis(LOW_CONFIDENCE_TARGET['name'], hp_low, max_qubits=10)  # Can handle full 10 residues
    print("Results low done!")
    results_low_mj = run_vqe_analysis_mj(LOW_CONFIDENCE_TARGET['name'], LOW_CONFIDENCE_TARGET['sequence'] , energy_matrix, symbols, max_qubits=10)
    print("Results low - mj done!")
    
    # Summary
    print("\n" + "="*80)
    print("📊 QUANTUM ANALYSIS SUMMARY")
    print("="*80)
    print(f"\n🧬 {results_high['name']}")
    print(f"   Ground State Energy: {results_high['energy']:.4f}")
    print(f"   Stability Score: {results_high['stability_score']:.4f}")
    print(f"   Classification: {results_high['stability']}")
    print(f"   VQE Runtime: {results_high['time']:.2f}s")
    
    print(f"\n🧬 {results_low['name']}")
    print(f"   Ground State Energy: {results_low['energy']:.4f}")
    print(f"   Stability Score: {results_low['stability_score']:.4f}")
    print(f"   Classification: {results_low['stability']}")
    print(f"   VQE Runtime: {results_low['time']:.2f}s")
    
    print("\n" + "="*80)
    print("🎯 BIOLOGICAL INTERPRETATION")
    print("="*80)
    
    #if results_high['stability_score'] > results_low['stability_score']:
    #    print("\n✅ EXPECTED RESULT CONFIRMED!")
    #    print(f"   The structured RBD core shows HIGHER stability ({results_high['stability_score']:.2f})")
    #    print(f"   The flexible loop shows LOWER stability ({results_low['stability_score']:.2f})")
    #    print("\n   This quantum simulation confirms that:")
    #    print("   • The RBD core has favorable hydrophobic interactions")
    #    print("   • The flexible loop is less energetically stable")
    #    print("   • Structural predictions align with quantum energy landscape")
    #else:
    #    print("\n⚠️  UNEXPECTED RESULT")
    #    print("   Further analysis needed - could indicate:")
    #    print("   • Limitations of simplified HP model")
    #    print("   • Need for longer VQE optimization")
    #    print("   • Importance of longer-range interactions")
    #
    #print("\n" + "="*80)
    #print("📈 Generating visualizations...")
    #plot_comparison(results_high, results_low)

    if results_high_mj['stability_score'] > results_low_mj['stability_score']:
        print("\n✅ EXPECTED RESULT CONFIRMED!")
        print(f"   The structured RBD core shows HIGHER stability ({results_high_mj['stability_score']:.2f})")
        print(f"   The flexible loop shows LOWER stability ({results_low_mj['stability_score']:.2f})")
        print("\n   This quantum simulation confirms that:")
        print("   • The RBD core has favorable hydrophobic interactions")
        print("   • The flexible loop is less energetically stable")
        print("   • Structural predictions align with quantum energy landscape")
    else:
        print("\n⚠️  UNEXPECTED RESULT")
        print("   Further analysis needed - could indicate:")
        print("   • Limitations of simplified HP model")
        print("   • Need for longer VQE optimization")
        print("   • Importance of longer-range interactions")
    
    print("\n" + "="*80)
    print("📈 Generating visualizations...")
    plot_comparison(results_high_mj, results_low_mj)
    
    print("\n✅ Analysis complete!")
    print("="*80)

if __name__ == "__main__":
    main()