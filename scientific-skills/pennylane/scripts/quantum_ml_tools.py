#!/usr/bin/env python3
"""
PennyLane Tools - Quantum machine learning and circuit training utilities
Provides quantum circuit construction, training, and simulation tools
"""

import argparse
import json
import numpy as np
from typing import List, Callable, Optional


def create_variational_circuit(n_qubits: int, n_layers: int, dev_name: str = "default.qubit"):
    """Create a variational quantum circuit."""
    try:
        import pennylane as qml
    except ImportError:
        print("Error: PennyLane not installed. Run: uv pip install pennylane")
        return None, None
    
    dev = qml.device(dev_name, wires=n_qubits)
    
    @qml.qnode(dev)
    def circuit(params, x=None):
        # Data encoding
        if x is not None:
            qml.AngleEmbedding(x, wires=range(n_qubits))
        
        # Variational layers
        for layer in range(n_layers):
            for i in range(n_qubits):
                qml.RX(params[layer, i, 0], wires=i)
                qml.RY(params[layer, i, 1], wires=i)
                qml.RZ(params[layer, i, 2], wires=i)
            
            # Entanglement
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i + 1])
        
        return qml.expval(qml.PauliZ(0))
    
    return circuit, dev


def train_vqe(hamiltonian, n_qubits: int, n_layers: int = 2, steps: int = 100):
    """Train Variational Quantum Eigensolver."""
    try:
        import pennylane as qml
    except ImportError:
        print("Error: PennyLane not installed")
        return None
    
    dev = qml.device("default.qubit", wires=n_qubits)
    
    @qml.qnode(dev)
    def vqe_circuit(params):
        for layer in range(n_layers):
            for i in range(n_qubits):
                qml.RY(params[layer, i], wires=i)
            for i in range(n_qubits - 1):
                qml.CNOT(wires=[i, i + 1])
        return qml.expval(hamiltonian)
    
    # Initialize parameters
    params = np.random.uniform(0, 2 * np.pi, (n_layers, n_qubits))
    opt = qml.GradientDescentOptimizer(stepsize=0.1)
    
    energies = []
    for i in range(steps):
        params, energy = opt.step_and_cost(vqe_circuit, params)
        energies.append(float(energy))
        if i % 10 == 0:
            print(f"Step {i}: Energy = {energy:.6f}")
    
    return {
        "final_energy": float(energy),
        "energies": energies,
        "optimal_params": params.tolist()
    }


def create_quantum_classifier(n_features: int, n_classes: int = 2, n_layers: int = 3):
    """Create a quantum classifier."""
    try:
        import pennylane as qml
    except ImportError:
        print("Error: PennyLane not installed")
        return None, None
    
    n_qubits = int(np.ceil(np.log2(n_classes))) if n_classes > 2 else 1
    n_qubits = max(n_qubits, n_features)
    
    dev = qml.device("default.qubit", wires=n_qubits)
    
    @qml.qnode(dev)
    def classifier(x, weights):
        # Angle embedding
        qml.AngleEmbedding(x, wires=range(n_qubits))
        
        # Variational layers
        qml.StronglyEntanglingLayers(weights, wires=range(n_qubits))
        
        # Measurement
        if n_classes == 2:
            return qml.expval(qml.PauliZ(0))
        else:
            return [qml.expval(qml.PauliZ(i)) for i in range(n_qubits)]
    
    return classifier, n_qubits


def visualize_circuit(circuit_func, params, filename: str = "circuit.png"):
    """Visualize a quantum circuit."""
    try:
        import pennylane as qml
        import matplotlib.pyplot as plt
        
        fig, ax = qml.draw_mpl(circuit_func)(params)
        plt.savefig(filename)
        print(f"Circuit saved to {filename}")
        return True
    except Exception as e:
        print(f"Error visualizing circuit: {e}")
        return False


def benchmark_devices(n_qubits: int, n_layers: int = 3, n_runs: int = 10):
    """Benchmark different PennyLane devices."""
    try:
        import pennylane as qml
        import time
    except ImportError:
        print("Error: PennyLane not installed")
        return {}
    
    devices = [
        ("default.qubit", {"wires": n_qubits}),
        ("lightning.qubit", {"wires": n_qubits}),
    ]
    
    results = {}
    
    for dev_name, dev_kwargs in devices:
        try:
            dev = qml.device(dev_name, **dev_kwargs)
            
            @qml.qnode(dev)
            def benchmark_circuit(params):
                for layer in range(n_layers):
                    for i in range(n_qubits):
                        qml.RX(params[layer, i, 0], wires=i)
                        qml.RY(params[layer, i, 1], wires=i)
                    for i in range(n_qubits - 1):
                        qml.CNOT(wires=[i, i + 1])
                return qml.expval(qml.PauliZ(0))
            
            params = np.random.random((n_layers, n_qubits, 2))
            
            # Warmup
            benchmark_circuit(params)
            
            # Benchmark
            start = time.time()
            for _ in range(n_runs):
                benchmark_circuit(params)
            elapsed = time.time() - start
            
            results[dev_name] = {
                "total_time": elapsed,
                "per_run": elapsed / n_runs
            }
            print(f"{dev_name}: {elapsed:.4f}s total, {elapsed/n_runs:.4f}s per run")
            
        except Exception as e:
            print(f"{dev_name} failed: {e}")
            results[dev_name] = {"error": str(e)}
    
    return results


def export_circuit_qasm(circuit_func, params, filename: str = "circuit.qasm"):
    """Export circuit to OpenQASM format."""
    try:
        import pennylane as qml
        
        qasm_str = qml.qasm(circuit_func)(params)
        
        with open(filename, 'w') as f:
            f.write(qasm_str)
        
        print(f"QASM exported to {filename}")
        return True
    except Exception as e:
        print(f"Error exporting QASM: {e}")
        return False


def create_h2_hamiltonian(bond_length: float = 0.74):
    """Create H2 molecule Hamiltonian."""
    try:
        import pennylane as qml
        from pennylane import qchem
        
        symbols = ["H", "H"]
        coordinates = np.array([0.0, 0.0, 0.0, 0.0, 0.0, bond_length])
        
        H, qubits = qchem.molecular_hamiltonian(symbols, coordinates)
        
        print(f"H2 Hamiltonian created ({qubits} qubits)")
        print(f"Terms: {len(H.ops)}")
        
        return H, qubits
    except Exception as e:
        print(f"Error creating Hamiltonian: {e}")
        return None, 0


def main():
    parser = argparse.ArgumentParser(description="PennyLane Quantum ML Tools")
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # VQE command
    vqe_parser = subparsers.add_parser("vqe", help="Run VQE for H2 molecule")
    vqe_parser.add_argument("-b", "--bond-length", type=float, default=0.74)
    vqe_parser.add_argument("-l", "--layers", type=int, default=2)
    vqe_parser.add_argument("-s", "--steps", type=int, default=100)
    vqe_parser.add_argument("-o", "--output", help="Output JSON file")
    
    # Benchmark command
    bench_parser = subparsers.add_parser("benchmark", help="Benchmark devices")
    bench_parser.add_argument("-q", "--qubits", type=int, default=4)
    bench_parser.add_argument("-l", "--layers", type=int, default=3)
    bench_parser.add_argument("-r", "--runs", type=int, default=10)
    
    # Hamiltonian command
    ham_parser = subparsers.add_parser("hamiltonian", help="Create molecular Hamiltonian")
    ham_parser.add_argument("-b", "--bond-length", type=float, default=0.74)
    
    # Circuit spec command
    spec_parser = subparsers.add_parser("specs", help="Get circuit specifications")
    spec_parser.add_argument("-q", "--qubits", type=int, default=4)
    spec_parser.add_argument("-l", "--layers", type=int, default=3)
    
    args = parser.parse_args()
    
    if args.command == "vqe":
        H, n_qubits = create_h2_hamiltonian(args.bond_length)
        if H is not None:
            results = train_vqe(H, n_qubits, args.layers, args.steps)
            if args.output:
                with open(args.output, 'w') as f:
                    json.dump(results, f, indent=2)
                print(f"Results saved to {args.output}")
    
    elif args.command == "benchmark":
        results = benchmark_devices(args.qubits, args.layers, args.runs)
        print(json.dumps(results, indent=2))
    
    elif args.command == "hamiltonian":
        H, n_qubits = create_h2_hamiltonian(args.bond_length)
        if H is not None:
            print(f"\nHamiltonian: {H}")
    
    elif args.command == "specs":
        try:
            import pennylane as qml
            
            dev = qml.device("default.qubit", wires=args.qubits)
            
            @qml.qnode(dev)
            def circuit(params):
                for layer in range(args.layers):
                    for i in range(args.qubits):
                        qml.RX(params[layer, i, 0], wires=i)
                        qml.RY(params[layer, i, 1], wires=i)
                    for i in range(args.qubits - 1):
                        qml.CNOT(wires=[i, i + 1])
                return qml.expval(qml.PauliZ(0))
            
            params = np.random.random((args.layers, args.qubits, 2))
            specs = qml.specs(circuit)(params)
            
            print("Circuit Specifications:")
            print(json.dumps(specs, indent=2, default=str))
        except ImportError:
            print("PennyLane not installed")
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
