"""
FluidSim Example: Parametric Study
Run multiple simulations with different parameters.
"""
from fluidsim.solvers.ns2d.solver import Simul
from math import pi


def parametric_viscosity_study():
    """
    Run simulations with different viscosities.
    """
    viscosities = [1e-3, 5e-4, 1e-4]
    
    for nu in viscosities:
        print(f"\nRunning simulation with nu = {nu}")
        
        params = Simul.create_default_params()
        params.oper.nx = params.oper.ny = 256
        params.oper.Lx = params.oper.Ly = 2 * pi
        params.nu_2 = nu
        params.time_stepping.t_end = 5.0
        params.init_fields.type = "noise"
        
        # Organize output by parameter
        params.output.sub_directory = f"nu_{nu}"
        
        sim = Simul(params)
        sim.time_stepping.start()
        
        print(f"Completed: nu = {nu}")


def resolution_convergence_study():
    """
    Study convergence with different resolutions.
    """
    resolutions = [64, 128, 256]
    
    for nx in resolutions:
        print(f"\nRunning with resolution {nx}x{nx}")
        
        params = Simul.create_default_params()
        params.oper.nx = params.oper.ny = nx
        params.oper.Lx = params.oper.Ly = 2 * pi
        params.nu_2 = 1e-3
        params.time_stepping.t_end = 2.0
        params.init_fields.type = "noise"
        
        params.output.sub_directory = f"nx_{nx}"
        
        sim = Simul(params)
        sim.time_stepping.start()
        
        print(f"Completed: nx = {nx}")


if __name__ == "__main__":
    print("FluidSim Parametric Study Example")
    # parametric_viscosity_study()
    # resolution_convergence_study()
