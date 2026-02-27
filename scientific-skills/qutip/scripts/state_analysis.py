#!/usr/bin/env python3
"""
QuTiP State Analysis Tools
Utilities for analyzing quantum states and operators.
"""

import argparse
import json
import numpy as np


def analyze_state(state_vector: list = None, density_matrix: list = None) -> dict:
    """
    Analyze a quantum state and compute various properties.
    
    Args:
        state_vector: State vector (ket) as list
        density_matrix: Density matrix as 2D list
        
    Returns:
        Dictionary of state properties
    """
    try:
        import qutip
    except ImportError:
        raise ImportError("qutip is required. Install with: uv pip install qutip")
    
    # Create state object
    if state_vector is not None:
        psi = qutip.Qobj(np.array(state_vector))
        rho = psi.proj()  # Create density matrix
    elif density_matrix is not None:
        rho = qutip.Qobj(np.array(density_matrix))
        psi = None
    else:
        raise ValueError("Either state_vector or density_matrix must be provided")
    
    # Calculate properties
    analysis = {
        'is_pure': qutip.entropy_vn(rho) < 1e-10,
        'von_neumann_entropy': qutip.entropy_vn(rho),
    }
    
    # For two-level systems
    if rho.shape[0] == 2:
        # Bloch vector
        sx = qutip.expect(qutip.sigmax(), rho)
        sy = qutip.expect(qutip.sigmay(), rho)
        sz = qutip.expect(qutip.sigmaz(), rho)
        analysis['bloch_vector'] = [sx, sy, sz]
        analysis['bloch_magnitude'] = np.sqrt(sx**2 + sy**2 + sz**2)
    
    # For two-qubit systems
    if rho.shape[0] == 4:
        analysis['concurrence'] = qutip.concurrence(rho)
        analysis['purity'] = (rho * rho).tr().real
    
    return analysis


def compute_wigner_function(state_vector: list, xvec: list = None, 
                             num_points: int = 200) -> dict:
    """
    Compute Wigner function for a quantum state.
    
    Args:
        state_vector: State vector
        xvec: X coordinates (auto-generated if None)
        num_points: Number of points per dimension
        
    Returns:
        Wigner function data
    """
    try:
        import qutip
    except ImportError:
        raise ImportError("qutip is required")
    
    psi = qutip.Qobj(np.array(state_vector))
    
    # Generate coordinate grid
    if xvec is None:
        lim = 5.0
        xvec = np.linspace(-lim, lim, num_points)
    else:
        xvec = np.array(xvec)
    
    # Compute Wigner function
    W = qutip.wigner(psi, xvec, xvec)
    
    return {
        'x': xvec.tolist(),
        'y': xvec.tolist(),
        'W': W.tolist(),
    }


def compare_states(state1: list, state2: list) -> dict:
    """
    Compare two quantum states.
    
    Args:
        state1: First state vector
        state2: Second state vector
        
    Returns:
        Comparison metrics
    """
    try:
        import qutip
    except ImportError:
        raise ImportError("qutip is required")
    
    psi1 = qutip.Qobj(np.array(state1))
    psi2 = qutip.Qobj(np.array(state2))
    
    rho1 = psi1.proj()
    rho2 = psi2.proj()
    
    return {
        'fidelity': qutip.fidelity(psi1, psi2),
        'tracedist': qutip.tracedist(rho1, rho2),
        'overlap': np.abs((psi1.dag() * psi2).full()[0, 0])**2,
    }


def generate_coherent_states(alpha_values: list, N: int = 20) -> dict:
    """
    Generate coherent states for given alpha values.
    
    Args:
        alpha_values: List of complex alpha values
        N: Fock space dimension
        
    Returns:
        Dictionary of coherent states
    """
    try:
        import qutip
    except ImportError:
        raise ImportError("qutip is required")
    
    states = {}
    for alpha in alpha_values:
        alpha_complex = complex(alpha)
        psi = qutip.coherent(N, alpha_complex)
        states[str(alpha)] = {
            'vector': psi.full().flatten().tolist(),
            'mean_photon': abs(alpha_complex)**2,
            'alpha': alpha,
        }
    
    return states


def generate_cat_states(alpha: float, N: int = 20) -> dict:
    """
    Generate cat states (superposition of coherent states).
    
    Args:
        alpha: Coherent state amplitude
        N: Fock space dimension
        
    Returns:
        Cat state data
    """
    try:
        import qutip
    except ImportError:
        raise ImportError("qutip is required")
    
    # Even cat state: |α> + |-α>
    psi_even = qutip.coherent(N, alpha) + qutip.coherent(N, -alpha)
    psi_even = psi_even.unit()
    
    # Odd cat state: |α> - |-α>
    psi_odd = qutip.coherent(N, alpha) - qutip.coherent(N, -alpha)
    psi_odd = psi_odd.unit()
    
    return {
        'even_cat': {
            'vector': psi_even.full().flatten().tolist(),
            'name': f'(|α> + |-α>), α={alpha}'
        },
        'odd_cat': {
            'vector': psi_odd.full().flatten().tolist(),
            'name': f'(|α> - |-α>), α={alpha}'
        },
        'alpha': alpha,
    }


def main():
    parser = argparse.ArgumentParser(description='QuTiP State Analysis Tools')
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Analyze command
    analyze_parser = subparsers.add_parser('analyze', help='Analyze quantum state')
    analyze_parser.add_argument('-s', '--state', nargs='+', type=float,
                                help='State vector as space-separated values')
    analyze_parser.add_argument('--output', help='Output JSON file')
    
    # Wigner command
    wigner_parser = subparsers.add_parser('wigner', help='Compute Wigner function')
    wigner_parser.add_argument('-s', '--state', nargs='+', type=float,
                               help='State vector')
    wigner_parser.add_argument('-n', '--points', type=int, default=100,
                               help='Number of grid points')
    wigner_parser.add_argument('--output', help='Output JSON file')
    
    # Coherent command
    coh_parser = subparsers.add_parser('coherent', help='Generate coherent states')
    coh_parser.add_argument('-a', '--alpha', nargs='+', type=float,
                            help='Alpha values (real)')
    coh_parser.add_argument('-N', type=int, default=20, help='Fock dimension')
    coh_parser.add_argument('--output', help='Output JSON file')
    
    # Cat state command
    cat_parser = subparsers.add_parser('cat', help='Generate cat states')
    cat_parser.add_argument('-a', '--alpha', type=float, default=2.0,
                            help='Alpha value')
    cat_parser.add_argument('-N', type=int, default=20, help='Fock dimension')
    cat_parser.add_argument('--output', help='Output JSON file')
    
    args = parser.parse_args()
    
    if args.command == 'analyze':
        if not args.state:
            print("Error: --state required")
            return
        results = analyze_state(state_vector=args.state)
        print("\nState Analysis:")
        for key, value in results.items():
            print(f"  {key}: {value}")
            
    elif args.command == 'wigner':
        if not args.state:
            print("Error: --state required")
            return
        results = compute_wigner_function(args.state, num_points=args.points)
        print(f"\nWigner function computed:")
        print(f"  Grid size: {args.points}x{args.points}")
        
    elif args.command == 'coherent':
        if not args.alpha:
            print("Error: --alpha required")
            return
        alphas = [complex(a, 0) for a in args.alpha]
        results = generate_coherent_states(alphas, args.N)
        print(f"\nGenerated {len(results)} coherent states")
        for alpha, data in results.items():
            print(f"  α={alpha}: <n> = {data['mean_photon']:.2f}")
            
    elif args.command == 'cat':
        results = generate_cat_states(args.alpha, args.N)
        print(f"\nGenerated cat states with α={args.alpha}")
        for name, data in results.items():
            if isinstance(data, dict) and 'name' in data:
                print(f"  {name}: {data['name']}")
            
    else:
        parser.print_help()
        return
    
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to: {args.output}")


if __name__ == '__main__':
    main()
