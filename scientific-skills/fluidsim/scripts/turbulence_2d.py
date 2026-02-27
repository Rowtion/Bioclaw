"""
FluidSim Example: 2D Turbulence Simulation
Run a 2D Navier-Stokes turbulence simulation.
"""
from fluidsim.solvers.ns2d.solver import Simul
from math import pi


def run_2d_turbulence():
    """
    Run a basic 2D turbulence simulation.
    """
    # Create default parameters
    params = Simul.create_default_params()
    
    # Set resolution and domain
    params.oper.nx = params.oper.ny = 256
    params.oper.Lx = params.oper.Ly = 2 * pi
    
    # Physical parameters
    params.nu_2 = 1e-3  # Viscosity
    
    # Time stepping
    params.time_stepping.t_end = 10.0
    params.time_stepping.USE_CFL = True
    params.time_stepping.CFL = 0.5
    
    # Initial conditions
    params.init_fields.type = "noise"
    
    # Output settings
    params.output.periods_save.phys_fields = 1.0
    params.output.periods_save.spectra = 0.5
    params.output.periods_save.spatial_means = 0.1
    
    # Create simulation
    sim = Simul(params)
    
    # Run simulation
    print("Starting 2D turbulence simulation...")
    sim.time_stepping.start()
    
    return sim


def analyze_turbulence_output(sim):
    """
    Analyze simulation output.
    
    Args:
        sim: Simul object after simulation
    """
    # Plot physical fields
    sim.output.phys_fields.plot("vorticity")
    sim.output.phys_fields.plot("vx")
    
    # Plot spectra
    sim.output.spectra.plot1d()
    
    # Plot spatial means
    sim.output.spatial_means.plot()


if __name__ == "__main__":
    print("FluidSim 2D Turbulence Example")
    print("Note: This requires FluidSim installation and may take time to run.")
    # sim = run_2d_turbulence()
    # analyze_turbulence_output(sim)
