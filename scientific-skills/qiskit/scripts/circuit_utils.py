#!/usr/bin/env python3
"""
Qiskit Quantum Circuit Utilities
Utilities for building, simulating, and optimizing quantum circuits.
"""

import argparse
import json
import os
from pathlib import Path


def create_bell_circuit(num_qubits: int = 2, name: str = "bell_state") -> str:
    """
    Create a Bell state circuit.
    
    Args:
        num_qubits: Number of qubits (must be even)
        name: Circuit name
        
    Returns:
        Qiskit QuantumCircuit JSON representation
    """
    try:
        from qiskit import QuantumCircuit
    except ImportError:
        raise ImportError("qiskit is required. Install with: uv pip install qiskit")
    
    if num_qubits % 2 != 0:
        raise ValueError("Number of qubits must be even")
    
    qc = QuantumCircuit(num_qubits, num_qubits, name=name)
    
    # Create Bell pairs
    for i in range(0, num_qubits, 2):
        qc.h(i)
        qc.cx(i, i + 1)
    
    # Measure all
    qc.measure_all()
    
    # Return circuit as dictionary
    return {
        'name': name,
        'num_qubits': num_qubits,
        'circuit_qasm': qc.qasm(),
        'depth': qc.depth(),
        'gate_count': sum(qc.count_ops().values())
    }


def create_ghz_circuit(num_qubits: int, name: str = "ghz_state") -> dict:
    """
    Create a GHZ (Greenberger-Horne-Zeilinger) state circuit.
    
    Args:
        num_qubits: Number of qubits
        name: Circuit name
        
    Returns:
        Circuit information dictionary
    """
    try:
        from qiskit import QuantumCircuit
    except ImportError:
        raise ImportError("qiskit is required")
    
    qc = QuantumCircuit(num_qubits, num_qubits, name=name)
    
    # Create GHZ state: |0...0> + |1...1>
    qc.h(0)
    for i in range(num_qubits - 1):
        qc.cx(i, i + 1)
    
    qc.measure_all()
    
    return {
        'name': name,
        'num_qubits': num_qubits,
        'circuit_qasm': qc.qasm(),
        'depth': qc.depth(),
        'gate_count': sum(qc.count_ops().values())
    }


def create_qft_circuit(num_qubits: int, inverse: bool = False, name: str = "qft") -> dict:
    """
    Create a Quantum Fourier Transform circuit.
    
    Args:
        num_qubits: Number of qubits
        inverse: Create inverse QFT
        name: Circuit name
        
    Returns:
        Circuit information dictionary
    """
    try:
        from qiskit import QuantumCircuit
        from qiskit.circuit.library import QFT
    except ImportError:
        raise ImportError("qiskit is required")
    
    qc = QuantumCircuit(num_qubits, num_qubits, name=name)
    
    # Add QFT
    qft = QFT(num_qubits, inverse=inverse)
    qc.compose(qft, inplace=True)
    
    qc.measure_all()
    
    return {
        'name': name,
        'num_qubits': num_qubits,
        'inverse': inverse,
        'circuit_qasm': qc.qasm(),
        'depth': qc.depth(),
        'gate_count': sum(qc.count_ops().values())
    }


def simulate_circuit(circuit_qasm: str, shots: int = 1024) -> dict:
    """
    Simulate a quantum circuit using Qiskit's statevector simulator.
    
    Args:
        circuit_qasm: Circuit in QASM format
        shots: Number of shots to run
        
    Returns:
        Simulation results
    """
    try:
        from qiskit import QuantumCircuit
        from qiskit.primitives import StatevectorSampler
    except ImportError:
        raise ImportError("qiskit is required")
    
    # Load circuit from QASM
    qc = QuantumCircuit.from_qasm_str(circuit_qasm)
    
    # Run simulation
    sampler = StatevectorSampler()
    result = sampler.run([qc], shots=shots).result()
    
    # Extract counts
    counts = result[0].data.meas.get_counts()
    
    return {
        'shots': shots,
        'counts': counts,
        'total_counts': sum(counts.values()),
        'num_results': len(counts)
    }


def transpile_circuit(circuit_qasm: str, backend_name: str = None, 
                      optimization_level: int = 3) -> dict:
    """
    Transpile a quantum circuit for a specific backend.
    
    Args:
        circuit_qasm: Circuit in QASM format
        backend_name: Target backend (None for generic)
        optimization_level: Optimization level (0-3)
        
    Returns:
        Transpiled circuit information
    """
    try:
        from qiskit import QuantumCircuit, transpile
        from qiskit_ibm_runtime import QiskitRuntimeService
    except ImportError:
        raise ImportError("qiskit and qiskit-ibm-runtime are required")
    
    # Load circuit
    qc = QuantumCircuit.from_qasm_str(circuit_qasm)
    
    # Get backend if specified
    backend = None
    if backend_name:
        service = QiskitRuntimeService()
        backend = service.backend(backend_name)
    
    # Transpile
    qc_transpiled = transpile(qc, backend=backend, optimization_level=optimization_level)
    
    return {
        'original_depth': qc.depth(),
        'transpiled_depth': qc_transpiled.depth(),
        'original_gates': qc.count_ops(),
        'transpiled_gates': qc_transpiled.count_ops(),
        'circuit_qasm': qc_transpiled.qasm(),
        'backend': backend_name if backend_name else 'generic'
    }


def calculate_circuit_metrics(circuit_qasm: str) -> dict:
    """
    Calculate various metrics for a quantum circuit.
    
    Args:
        circuit_qasm: Circuit in QASM format
        
    Returns:
        Dictionary of metrics
    """
    try:
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import Statevector, Operator
    except ImportError:
        raise ImportError("qiskit is required")
    
    qc = QuantumCircuit.from_qasm_str(circuit_qasm)
    
    # Remove measurements for unitary analysis
    qc_no_meas = qc.remove_final_measurements(inplace=False)
    
    metrics = {
        'num_qubits': qc.num_qubits,
        'num_clbits': qc.num_clbits,
        'depth': qc.depth(),
        'total_gates': sum(qc.count_ops().values()),
        'gate_counts': dict(qc.count_ops()),
        'width': qc.width(),
    }
    
    # Try to get unitary (may fail for large circuits)
    try:
        if qc_no_meas.num_qubits <= 10:
            unitary = Operator(qc_no_meas)
            metrics['is_unitary'] = True
    except:
        metrics['is_unitary'] = False
    
    return metrics


def main():
    parser = argparse.ArgumentParser(description='Qiskit Quantum Circuit Utilities')
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Bell command
    bell_parser = subparsers.add_parser('bell', help='Create Bell state circuit')
    bell_parser.add_argument('-n', '--num-qubits', type=int, default=2,
                             help='Number of qubits')
    bell_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # GHZ command
    ghz_parser = subparsers.add_parser('ghz', help='Create GHZ state circuit')
    ghz_parser.add_argument('-n', '--num-qubits', type=int, required=True,
                            help='Number of qubits')
    ghz_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # QFT command
    qft_parser = subparsers.add_parser('qft', help='Create QFT circuit')
    qft_parser.add_argument('-n', '--num-qubits', type=int, required=True,
                            help='Number of qubits')
    qft_parser.add_argument('--inverse', action='store_true',
                            help='Create inverse QFT')
    qft_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # Simulate command
    sim_parser = subparsers.add_parser('simulate', help='Simulate circuit')
    sim_parser.add_argument('input', help='Input QASM file')
    sim_parser.add_argument('-s', '--shots', type=int, default=1024,
                            help='Number of shots')
    sim_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # Transpile command
    transpile_parser = subparsers.add_parser('transpile', help='Transpile circuit')
    transpile_parser.add_argument('input', help='Input QASM file')
    transpile_parser.add_argument('-b', '--backend', help='Target backend')
    transpile_parser.add_argument('-l', '--level', type=int, default=3,
                                  choices=[0, 1, 2, 3],
                                  help='Optimization level')
    transpile_parser.add_argument('-o', '--output', help='Output QASM file')
    
    # Metrics command
    metrics_parser = subparsers.add_parser('metrics', help='Calculate circuit metrics')
    metrics_parser.add_argument('input', help='Input QASM file')
    
    args = parser.parse_args()
    
    if args.command == 'bell':
        result = create_bell_circuit(args.num_qubits)
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(result, f, indent=2)
            print(f"Circuit saved to: {args.output}")
        print(f"\nBell Circuit ({args.num_qubits} qubits):")
        print(f"  Depth: {result['depth']}")
        print(f"  Gates: {result['gate_count']}")
        
    elif args.command == 'ghz':
        result = create_ghz_circuit(args.num_qubits)
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(result, f, indent=2)
        print(f"\nGHZ Circuit ({args.num_qubits} qubits):")
        print(f"  Depth: {result['depth']}")
        print(f"  Gates: {result['gate_count']}")
        
    elif args.command == 'qft':
        result = create_qft_circuit(args.num_qubits, args.inverse)
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(result, f, indent=2)
        print(f"\nQFT Circuit ({args.num_qubits} qubits):")
        print(f"  Depth: {result['depth']}")
        print(f"  Gates: {result['gate_count']}")
        
    elif args.command == 'simulate':
        with open(args.input) as f:
            qasm = f.read()
        result = simulate_circuit(qasm, args.shots)
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(result, f, indent=2)
        print(f"\nSimulation Results ({args.shots} shots):")
        for bitstring, count in sorted(result['counts'].items()):
            print(f"  {bitstring}: {count}")
            
    elif args.command == 'transpile':
        with open(args.input) as f:
            qasm = f.read()
        result = transpile_circuit(qasm, args.backend, args.level)
        if args.output:
            with open(args.output, 'w') as f:
                f.write(result['circuit_qasm'])
        print(f"\nTranspilation Results:")
        print(f"  Original depth: {result['original_depth']}")
        print(f"  Transpiled depth: {result['transpiled_depth']}")
        
    elif args.command == 'metrics':
        with open(args.input) as f:
            qasm = f.read()
        metrics = calculate_circuit_metrics(qasm)
        print(f"\nCircuit Metrics:")
        for key, value in metrics.items():
            print(f"  {key}: {value}")
            
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
