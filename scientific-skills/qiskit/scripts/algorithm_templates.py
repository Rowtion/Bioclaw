#!/usr/bin/env python3
"""
Qiskit Algorithm Templates
Pre-built quantum algorithm implementations.
"""

import argparse
import json
import numpy as np


def vqe_template(num_qubits: int, num_layers: int = 2) -> dict:
    """
    Create a Variational Quantum Eigensolver (VQE) ansatz template.
    
    Args:
        num_qubits: Number of qubits
        num_layers: Number of layers in the ansatz
        
    Returns:
        Circuit template information
    """
    try:
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import EfficientSU2
    except ImportError:
        raise ImportError("qiskit is required")
    
    # Create ansatz
    ansatz = EfficientSU2(num_qubits, reps=num_layers)
    
    return {
        'algorithm': 'VQE',
        'num_qubits': num_qubits,
        'num_parameters': len(ansatz.parameters),
        'circuit_qasm': ansatz.qasm(),
        'parameters': [str(p) for p in ansatz.parameters],
        'description': 'EfficientSU2 ansatz for VQE'
    }


def qaoa_template(num_nodes: int, num_layers: int = 1) -> dict:
    """
    Create a QAOA (Quantum Approximate Optimization Algorithm) template.
    
    Args:
        num_nodes: Number of nodes in the graph
        num_layers: Number of QAOA layers
        
    Returns:
        Circuit template information
    """
    try:
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import QAOAAnsatz
        from qiskit.quantum_info import SparsePauliOp
    except ImportError:
        raise ImportError("qiskit is required")
    
    # Create a simple MaxCut Hamiltonian
    edges = [(i, (i + 1) % num_nodes) for i in range(num_nodes)]
    
    # Build cost Hamiltonian
    pauli_list = []
    for i, j in edges:
        pauli_list.append(('I' * i + 'Z' + 'I' * (j - i - 1) + 'Z' + 'I' * (num_nodes - j - 1), 0.5))
    
    cost_hamiltonian = SparsePauliOp.from_list(pauli_list)
    
    # Create QAOA ansatz
    qaoa = QAOAAnsatz(cost_hamiltonian, reps=num_layers)
    
    return {
        'algorithm': 'QAOA',
        'num_qubits': num_nodes,
        'num_parameters': len(qaoa.parameters),
        'circuit_qasm': qaoa.qasm(),
        'description': f'QAOA for MaxCut on {num_nodes}-node ring graph'
    }


def grover_template(num_qubits: int, num_solutions: int = 1) -> dict:
    """
    Create a Grover's algorithm template.
    
    Args:
        num_qubits: Number of qubits
        num_solutions: Number of solutions to search for
        
    Returns:
        Circuit template information
    """
    try:
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import GroverOperator
    except ImportError:
        raise ImportError("qiskit is required")
    
    # Create oracle (marks state |11...1> as solution)
    oracle = QuantumCircuit(num_qubits)
    oracle.h(num_qubits - 1)
    oracle.mcx(list(range(num_qubits - 1)), num_qubits - 1)
    oracle.h(num_qubits - 1)
    
    # Create Grover operator
    grover_op = GroverOperator(oracle)
    
    # Calculate optimal iterations
    N = 2 ** num_qubits
    optimal_iterations = int(np.round(np.pi / 4 * np.sqrt(N / num_solutions)))
    
    # Build full circuit
    qc = QuantumCircuit(num_qubits, num_qubits)
    qc.h(range(num_qubits))
    
    for _ in range(optimal_iterations):
        qc.compose(grover_op, inplace=True)
    
    qc.measure_all()
    
    return {
        'algorithm': "Grover's Search",
        'num_qubits': num_qubits,
        'num_iterations': optimal_iterations,
        'circuit_qasm': qc.qasm(),
        'description': f"Grover's search for {num_solutions} solution(s) in {N} items"
    }


def deutsch_jozsa_template(n: int) -> dict:
    """
    Create a Deutsch-Jozsa algorithm template.
    
    Args:
        n: Number of input qubits
        
    Returns:
        Circuit template information
    """
    try:
        from qiskit import QuantumCircuit
    except ImportError:
        raise ImportError("qiskit is required")
    
    qc = QuantumCircuit(n + 1, n)
    
    # Initialize
    qc.x(n)
    qc.h(range(n + 1))
    
    # Oracle (example: balanced function)
    for i in range(n):
        qc.cx(i, n)
    
    # Final Hadamards
    qc.h(range(n))
    
    # Measure
    qc.measure(range(n), range(n))
    
    return {
        'algorithm': 'Deutsch-Jozsa',
        'num_qubits': n + 1,
        'circuit_qasm': qc.qasm(),
        'description': f'Deutsch-Jozsa algorithm for {n}-bit function'
    }


def bernstein_vazirani_template(n: int, secret_string: str = None) -> dict:
    """
    Create a Bernstein-Vazirani algorithm template.
    
    Args:
        n: Number of input qubits
        secret_string: Secret binary string (auto-generated if None)
        
    Returns:
        Circuit template information
    """
    try:
        from qiskit import QuantumCircuit
    except ImportError:
        raise ImportError("qiskit is required")
    
    if secret_string is None:
        secret_string = ''.join(np.random.choice(['0', '1']) for _ in range(n))
    
    qc = QuantumCircuit(n + 1, n)
    
    # Initialize
    qc.x(n)
    qc.h(range(n + 1))
    
    # Oracle: encode secret string
    for i, bit in enumerate(reversed(secret_string)):
        if bit == '1':
            qc.cx(i, n)
    
    # Final Hadamards
    qc.h(range(n))
    
    # Measure
    qc.measure(range(n), range(n))
    
    return {
        'algorithm': 'Bernstein-Vazirani',
        'num_qubits': n + 1,
        'secret_string': secret_string,
        'circuit_qasm': qc.qasm(),
        'description': f'Bernstein-Vazirani algorithm to find secret string {secret_string}'
    }


def simon_template(n: int) -> dict:
    """
    Create a Simon's algorithm template.
    
    Args:
        n: Number of input qubits
        
    Returns:
        Circuit template information
    """
    try:
        from qiskit import QuantumCircuit
    except ImportError:
        raise ImportError("qiskit is required")
    
    qc = QuantumCircuit(2 * n, n)
    
    # Initialize first register
    qc.h(range(n))
    
    # Oracle (example: f(x) = x xor s)
    for i in range(n):
        qc.cx(i, i + n)
    
    # Measure second register (for simulation only)
    qc.measure(range(n, 2 * n), range(n))
    
    # Hadamard on first register
    qc.h(range(n))
    
    # Measure first register
    qc.measure(range(n), range(n))
    
    return {
        'algorithm': "Simon's",
        'num_qubits': 2 * n,
        'circuit_qasm': qc.qasm(),
        'description': f"Simon's algorithm for period finding in {n}-bit function"
    }


def main():
    parser = argparse.ArgumentParser(description='Qiskit Algorithm Templates')
    subparsers = parser.add_subparsers(dest='command', help='Algorithms')
    
    # VQE command
    vqe_parser = subparsers.add_parser('vqe', help='VQE template')
    vqe_parser.add_argument('-n', '--num-qubits', type=int, required=True,
                            help='Number of qubits')
    vqe_parser.add_argument('-l', '--layers', type=int, default=2,
                            help='Number of layers')
    vqe_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # QAOA command
    qaoa_parser = subparsers.add_parser('qaoa', help='QAOA template')
    qaoa_parser.add_argument('-n', '--num-nodes', type=int, required=True,
                             help='Number of nodes')
    qaoa_parser.add_argument('-l', '--layers', type=int, default=1,
                             help='Number of layers')
    qaoa_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # Grover command
    grover_parser = subparsers.add_parser('grover', help="Grover's search template")
    grover_parser.add_argument('-n', '--num-qubits', type=int, required=True,
                               help='Number of qubits')
    grover_parser.add_argument('-s', '--num-solutions', type=int, default=1,
                               help='Number of solutions')
    grover_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # Deutsch-Jozsa command
    dj_parser = subparsers.add_parser('deutsch-jozsa', help='Deutsch-Jozsa template')
    dj_parser.add_argument('-n', type=int, required=True, help='Number of input qubits')
    dj_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # Bernstein-Vazirani command
    bv_parser = subparsers.add_parser('bernstein-vazirani', help='Bernstein-Vazirani template')
    bv_parser.add_argument('-n', type=int, required=True, help='Number of input qubits')
    bv_parser.add_argument('-s', '--secret', help='Secret binary string')
    bv_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # Simon command
    simon_parser = subparsers.add_parser('simon', help="Simon's algorithm template")
    simon_parser.add_argument('-n', type=int, required=True, help='Number of input qubits')
    simon_parser.add_argument('-o', '--output', help='Output JSON file')
    
    args = parser.parse_args()
    
    if args.command == 'vqe':
        result = vqe_template(args.num_qubits, args.layers)
        
    elif args.command == 'qaoa':
        result = qaoa_template(args.num_nodes, args.layers)
        
    elif args.command == 'grover':
        result = grover_template(args.num_qubits, args.num_solutions)
        
    elif args.command == 'deutsch-jozsa':
        result = deutsch_jozsa_template(args.n)
        
    elif args.command == 'bernstein-vazirani':
        result = bernstein_vazirani_template(args.n, args.secret)
        
    elif args.command == 'simon':
        result = simon_template(args.n)
        
    else:
        parser.print_help()
        return
    
    # Output results
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"Template saved to: {args.output}")
    
    print(f"\n{result['algorithm']} Template:")
    print(f"  Qubits: {result['num_qubits']}")
    if 'num_parameters' in result:
        print(f"  Parameters: {result['num_parameters']}")
    print(f"  Description: {result['description']}")


if __name__ == '__main__':
    main()
