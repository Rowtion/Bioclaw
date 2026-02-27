#!/usr/bin/env python3
"""
QuTiP Quantum Dynamics Simulator
Utilities for simulating quantum system dynamics and open quantum systems.
"""

import argparse
import json
import numpy as np


def simulate_rabi_oscillations(omega: float, t_max: float, num_points: int = 1000) -> dict:
    """
    Simulate Rabi oscillations for a two-level system.
    
    Args:
        omega: Rabi frequency
        t_max: Maximum simulation time
        num_points: Number of time points
        
    Returns:
        Simulation results with time and population data
    """
    try:
        import qutip
    except ImportError:
        raise ImportError("qutip is required. Install with: uv pip install qutip")
    
    # Two-level system Hamiltonian
    H = omega * qutip.sigmax()
    
    # Initial state: |0>
    psi0 = qutip.basis(2, 0)
    
    # Time list
    tlist = np.linspace(0, t_max, num_points)
    
    # Expectation values: <σz> (population inversion)
    result = qutip.sesolve(H, psi0, tlist, e_ops=[qutip.sigmaz(), qutip.sigmax(), qutip.sigmay()])
    
    return {
        'system': 'Two-level Rabi oscillations',
        'rabi_frequency': omega,
        'time': tlist.tolist(),
        'expectation_sz': result.expect[0].tolist(),
        'expectation_sx': result.expect[1].tolist(),
        'expectation_sy': result.expect[2].tolist(),
    }


def simulate_damped_oscillator(omega: float, kappa: float, alpha0: float,
                               num_fock: int = 20, t_max: float = 50.0) -> dict:
    """
    Simulate a damped quantum harmonic oscillator.
    
    Args:
        omega: Oscillator frequency
        kappa: Decay rate
        alpha0: Initial coherent state amplitude
        num_fock: Fock space dimension
        t_max: Maximum simulation time
        
    Returns:
        Simulation results
    """
    try:
        import qutip
    except ImportError:
        raise ImportError("qutip is required")
    
    # Operators
    a = qutip.destroy(num_fock)
    H = omega * a.dag() * a
    
    # Collapse operators (damping)
    c_ops = [np.sqrt(kappa) * a]
    
    # Initial state: coherent state
    psi0 = qutip.coherent(num_fock, alpha0)
    
    # Time evolution
    tlist = np.linspace(0, t_max, 500)
    result = qutip.mesolve(H, psi0, tlist, c_ops, e_ops=[a.dag() * a, a])
    
    return {
        'system': 'Damped harmonic oscillator',
        'omega': omega,
        'kappa': kappa,
        'alpha0': alpha0,
        'time': tlist.tolist(),
        'photon_number': result.expect[0].tolist(),
        'coherent_amplitude': np.abs(result.expect[1]).tolist(),
    }


def simulate_jaynes_cummings(N: int, wc: float, wa: float, g: float,
                              kappa: float, gamma: float, t_max: float = 50.0) -> dict:
    """
    Simulate Jaynes-Cummings model (cavity QED).
    
    Args:
        N: Cavity Fock space dimension
        wc: Cavity frequency
        wa: Atom frequency
        g: Coupling strength
        kappa: Cavity decay rate
        gamma: Atomic decay rate
        t_max: Maximum simulation time
        
    Returns:
        Simulation results
    """
    try:
        import qutip
    except ImportError:
        raise ImportError("qutip is required")
    
    # Operators
    a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
    sm = qutip.tensor(qutip.qeye(N), qutip.sigmam())
    
    # Hamiltonian (RWA)
    H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())
    
    # Collapse operators
    c_ops = []
    if kappa > 0:
        c_ops.append(np.sqrt(kappa) * a)
    if gamma > 0:
        c_ops.append(np.sqrt(gamma) * sm)
    
    # Initial state: cavity in coherent state, atom in ground state
    psi0 = qutip.tensor(qutip.coherent(N, 2), qutip.basis(2, 0))
    
    # Observables
    n_cav = a.dag() * a
    n_atom = sm.dag() * sm
    
    # Evolve
    tlist = np.linspace(0, t_max, 500)
    result = qutip.mesolve(H, psi0, tlist, c_ops, e_ops=[n_cav, n_atom])
    
    return {
        'system': 'Jaynes-Cummings',
        'N': N,
        'wc': wc,
        'wa': wa,
        'g': g,
        'kappa': kappa,
        'gamma': gamma,
        'time': tlist.tolist(),
        'cavity_photons': result.expect[0].tolist(),
        'atom_excited': result.expect[1].tolist(),
    }


def simulate_entanglement_decay(gamma: float, t_max: float = 10.0) -> dict:
    """
    Simulate entanglement decay in a two-qubit system.
    
    Args:
        gamma: Dephasing rate
        t_max: Maximum simulation time
        
    Returns:
        Simulation results with concurrence
    """
    try:
        import qutip
    except ImportError:
        raise ImportError("qutip is required")
    
    # Initial Bell state
    psi0 = qutip.bell_state('00')
    
    # Local dephasing on each qubit
    c_ops = [
        np.sqrt(gamma) * qutip.tensor(qutip.sigmaz(), qutip.qeye(2)),
        np.sqrt(gamma) * qutip.tensor(qutip.qeye(2), qutip.sigmaz())
    ]
    
    # Time evolution (no Hamiltonian, just dissipation)
    tlist = np.linspace(0, t_max, 200)
    result = qutip.mesolve(qutip.qeye([2, 2]), psi0, tlist, c_ops)
    
    # Calculate concurrence for each state
    C_t = [qutip.concurrence(state.proj()) for state in result.states]
    
    return {
        'system': 'Two-qubit entanglement decay',
        'gamma': gamma,
        'time': tlist.tolist(),
        'concurrence': C_t,
    }


def simulate_spontaneous_emission(gamma: float, t_max: float = 10.0) -> dict:
    """
    Simulate spontaneous emission from excited state.
    
    Args:
        gamma: Decay rate
        t_max: Maximum simulation time
        
    Returns:
        Simulation results
    """
    try:
        import qutip
    except ImportError:
        raise ImportError("qutip is required")
    
    # Two-level system
    H = 0 * qutip.sigmaz()  # No coherent evolution
    
    # Initial state: excited
    psi0 = qutip.basis(2, 1)
    
    # Collapse operator
    c_ops = [np.sqrt(gamma) * qutip.sigmam()]
    
    # Observables
    e_ops = [qutip.ket2dm(qutip.basis(2, 1)), qutip.ket2dm(qutip.basis(2, 0))]
    
    # Evolve
    tlist = np.linspace(0, t_max, 200)
    result = qutip.mesolve(H, psi0, tlist, c_ops, e_ops=e_ops)
    
    return {
        'system': 'Spontaneous emission',
        'gamma': gamma,
        'time': tlist.tolist(),
        'excited_population': result.expect[0].tolist(),
        'ground_population': result.expect[1].tolist(),
    }


def save_results(results: dict, output_file: str):
    """Save simulation results to JSON file."""
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Results saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='QuTiP Quantum Dynamics Simulator')
    subparsers = parser.add_subparsers(dest='command', help='Simulations')
    
    # Rabi command
    rabi_parser = subparsers.add_parser('rabi', help='Rabi oscillations')
    rabi_parser.add_argument('-o', '--omega', type=float, default=1.0,
                             help='Rabi frequency')
    rabi_parser.add_argument('-t', '--tmax', type=float, default=10.0,
                             help='Max time')
    rabi_parser.add_argument('--output', help='Output JSON file')
    
    # Damped oscillator command
    damp_parser = subparsers.add_parser('damped', help='Damped oscillator')
    damp_parser.add_argument('--omega', type=float, default=1.0, help='Frequency')
    damp_parser.add_argument('--kappa', type=float, default=0.1, help='Decay rate')
    damp_parser.add_argument('--alpha', type=float, default=3.0, help='Initial amplitude')
    damp_parser.add_argument('--output', help='Output JSON file')
    
    # JC command
    jc_parser = subparsers.add_parser('jc', help='Jaynes-Cummings')
    jc_parser.add_argument('-N', type=int, default=10, help='Fock space dimension')
    jc_parser.add_argument('--wc', type=float, default=1.0, help='Cavity frequency')
    jc_parser.add_argument('--wa', type=float, default=1.0, help='Atom frequency')
    jc_parser.add_argument('-g', type=float, default=0.05, help='Coupling')
    jc_parser.add_argument('--kappa', type=float, default=0.1, help='Cavity decay')
    jc_parser.add_argument('--gamma', type=float, default=0.05, help='Atomic decay')
    jc_parser.add_argument('--output', help='Output JSON file')
    
    # Entanglement decay command
    ent_parser = subparsers.add_parser('entanglement', help='Entanglement decay')
    ent_parser.add_argument('-g', '--gamma', type=float, default=0.1, help='Dephasing rate')
    ent_parser.add_argument('--output', help='Output JSON file')
    
    # Spontaneous emission command
    se_parser = subparsers.add_parser('emission', help='Spontaneous emission')
    se_parser.add_argument('-g', '--gamma', type=float, default=1.0, help='Decay rate')
    se_parser.add_argument('--output', help='Output JSON file')
    
    args = parser.parse_args()
    
    if args.command == 'rabi':
        results = simulate_rabi_oscillations(args.omega, args.tmax)
        print(f"Rabi oscillations simulated:")
        print(f"  ω = {args.omega}")
        print(f"  Time points: {len(results['time'])}")
        
    elif args.command == 'damped':
        results = simulate_damped_oscillator(args.omega, args.kappa, args.alpha)
        print(f"Damped oscillator simulated:")
        print(f"  Initial photons: {args.alpha**2:.1f}")
        
    elif args.command == 'jc':
        results = simulate_jaynes_cummings(args.N, args.wc, args.wa, args.g,
                                           args.kappa, args.gamma)
        print(f"Jaynes-Cummings simulated:")
        print(f"  Cavity: N={args.N}, ωc={args.wc}")
        print(f"  Atom: ωa={args.wa}, g={args.g}")
        
    elif args.command == 'entanglement':
        results = simulate_entanglement_decay(args.gamma)
        print(f"Entanglement decay simulated:")
        print(f"  Initial concurrence: {results['concurrence'][0]:.3f}")
        print(f"  Final concurrence: {results['concurrence'][-1]:.3f}")
        
    elif args.command == 'emission':
        results = simulate_spontaneous_emission(args.gamma)
        print(f"Spontaneous emission simulated:")
        print(f"  Decay rate γ = {args.gamma}")
        
    else:
        parser.print_help()
        return
    
    if args.output:
        save_results(results, args.output)


if __name__ == '__main__':
    main()
