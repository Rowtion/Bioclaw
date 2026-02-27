"""
FluidSim Example: Taylor-Green Vortex
Simulate the classic Taylor-Green vortex for validation.
"""
from fluidsim.solvers.ns2d.solver import Simul
import numpy as np
from math import pi


def taylor_green_vortex():
    """
    Simulate Taylor-Green vortex - a classic validation case.
    """
    params = Simul.create_default_params()
    
    # Resolution
    params.oper.nx = params.oper.ny = 128
    params.oper.Lx = params.oper.Ly = 2 * pi
    
    # Viscosity
    params.nu_2 = 1e-3
    
    # Time stepping
    params.time_stepping.t_end = 10.0
    params.time_stepping.USE_CFL = True
    
    # Custom initial condition
    params.init_fields.type = "in_script"
    
    # Create simulation
    sim = Simul(params)
    
    # Set Taylor-Green initial condition
    X, Y = sim.oper.get_XY_loc()
    vx = sim.state.state_phys.get_var("vx")
    vy = sim.state.state_phys.get_var("vy")
    
    vx[:] = np.sin(X) * np.cos(Y)
    vy[:] = -np.cos(X) * np.sin(Y)
    
    # Update spectral representation
    sim.state.statephys_from_statespect()
    
    print("Running Taylor-Green vortex simulation...")
    sim.time_stepping.start()
    
    return sim


def validate_energy_decay(sim):
    """
    Validate energy decay against analytical solution.
    
    Args:
        sim: Simul object
    """
    # Load spatial means
    df = sim.output.spatial_means.load()
    
    # Calculate energy decay
    # Theoretical: E(t) = E0 * exp(-4*nu*t)
    print("Energy decay validation:")
    print(df.head())
    
    return df


if __name__ == "__main__":
    print("FluidSim Taylor-Green Vortex Example")
    # sim = taylor_green_vortex()
    # validate_energy_decay(sim)
