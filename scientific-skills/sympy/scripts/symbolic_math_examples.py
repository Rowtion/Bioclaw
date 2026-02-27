#!/usr/bin/env python3
"""
SymPy: Symbolic Mathematics Examples
=====================================
Exact symbolic computation for algebra, calculus, and more.
"""

import sympy as sp
from sympy import symbols, Symbol, simplify, expand, factor, diff, integrate, limit, oo
from sympy import sin, cos, tan, exp, log, sqrt, pi, E, I, Matrix, Rational
import numpy as np
import matplotlib.pyplot as plt

# Enable pretty printing
sp.init_printing()

# ==============================================================================
# Example 1: Symbolic Expressions and Simplification
# ==============================================================================

def example_expressions():
    """Working with symbolic expressions."""
    print("=" * 60)
    print("Example 1: Symbolic Expressions")
    print("=" * 60)
    
    # Define symbols
    x, y, z = symbols('x y z')
    
    # Create expressions
    expr1 = x**2 + 2*x + 1
    expr2 = (x + 1)**2
    expr3 = x**2 - y**2
    
    print(f"Expression 1: {expr1}")
    print(f"Expression 2: {expr2}")
    print(f"Expression 3: {expr3}")
    
    # Simplification
    print(f"\nSimplification:")
    print(f"  simplify({expr1} - {expr2}) = {simplify(expr1 - expr2)}")
    print(f"  factor({expr3}) = {factor(expr3)}")
    print(f"  expand({expr2}) = {expand(expr2)}")
    
    # More complex simplification
    expr4 = sin(x)**2 + cos(x)**2
    print(f"\n  simplify(sin(x)^2 + cos(x)^2) = {simplify(expr4)}")
    
    expr5 = (x**2 - 1) / (x - 1)
    print(f"  simplify((x^2-1)/(x-1)) = {simplify(expr5)}")
    
    # Trigonometric simplification
    expr6 = sin(2*x)
    print(f"\nTrigonometric:")
    print(f"  expand_trig(sin(2x)) = {sp.expand_trig(expr6)}")
    
    return expr1, expr2, expr3


# ==============================================================================
# Example 2: Equation Solving
# ==============================================================================

def example_solving():
    """Solving algebraic and transcendental equations."""
    print("\n" + "=" * 60)
    print("Example 2: Equation Solving")
    print("=" * 60)
    
    x, y = symbols('x y')
    
    # Quadratic equation
    eq1 = x**2 - 5*x + 6
    solutions1 = sp.solve(eq1, x)
    print(f"\nSolve x^2 - 5x + 6 = 0:")
    print(f"  Solutions: {solutions1}")
    
    # Verify solutions
    for sol in solutions1:
        result = eq1.subs(x, sol)
        print(f"  Substitute x={sol}: {simplify(result)}")
    
    # Cubic equation
    eq2 = x**3 - 6*x**2 + 11*x - 6
    solutions2 = sp.solve(eq2, x)
    print(f"\nSolve x^3 - 6x^2 + 11x - 6 = 0:")
    print(f"  Solutions: {solutions2}")
    
    # System of linear equations
    print(f"\nSystem of linear equations:")
    print(f"  x + y = 5")
    print(f"  x - y = 1")
    sol_system = sp.linsolve([x + y - 5, x - y - 1], x, y)
    print(f"  Solution: {sol_system}")
    
    # Nonlinear system
    print(f"\nNonlinear system:")
    print(f"  x^2 + y = 2")
    print(f"  x + y^2 = 3")
    sol_nonlinear = sp.nonlinsolve([x**2 + y - 2, x + y**2 - 3], x, y)
    print(f"  Solutions: {sol_nonlinear}")
    
    # Transcendental equation
    print(f"\nTranscendental equation:")
    eq3 = sp.Eq(exp(x), 10)
    sol_exp = sp.solve(eq3, x)
    print(f"  Solve e^x = 10: x = {sol_exp}")
    
    return solutions1


# ==============================================================================
# Example 3: Calculus - Derivatives
# ==============================================================================

def example_derivatives():
    """Symbolic differentiation."""
    print("\n" + "=" * 60)
    print("Example 3: Derivatives")
    print("=" * 60)
    
    x, y = symbols('x y')
    
    # Basic derivatives
    print(f"\nBasic derivatives:")
    print(f"  d/dx(x^3) = {diff(x**3, x)}")
    print(f"  d/dx(sin(x)) = {diff(sin(x), x)}")
    print(f"  d/dx(e^x) = {diff(exp(x), x)}")
    print(f"  d/dx(ln(x)) = {diff(log(x), x)}")
    
    # Higher-order derivatives
    print(f"\nHigher-order derivatives:")
    print(f"  d^2/dx^2(x^4) = {diff(x**4, x, 2)}")
    print(f"  d^3/dx^3(x^5) = {diff(x**5, x, 3)}")
    
    # Chain rule
    print(f"\nChain rule:")
    print(f"  d/dx(sin(x^2)) = {diff(sin(x**2), x)}")
    print(f"  d/dx(e^(3x)) = {diff(exp(3*x), x)}")
    
    # Product rule
    expr = x**2 * sin(x)
    print(f"\nProduct rule:")
    print(f"  d/dx(x^2 * sin(x)) = {diff(expr, x)}")
    
    # Quotient rule
    expr2 = sin(x) / x
    print(f"\nQuotient rule:")
    print(f"  d/dx(sin(x)/x) = {diff(expr2, x)}")
    
    # Partial derivatives
    f = x**2 * y + x * y**3
    print(f"\nPartial derivatives of f = x^2*y + x*y^3:")
    print(f"  ∂f/∂x = {diff(f, x)}")
    print(f"  ∂f/∂y = {diff(f, y)}")
    print(f"  ∂²f/∂x∂y = {diff(f, x, y)}")
    
    return diff(x**3, x)


# ==============================================================================
# Example 4: Calculus - Integrals
# ==============================================================================

def example_integrals():
    """Symbolic integration."""
    print("\n" + "=" * 60)
    print("Example 4: Integrals")
    print("=" * 60)
    
    x = symbols('x')
    
    # Indefinite integrals
    print(f"\nIndefinite integrals:")
    print(f"  ∫ x^2 dx = {integrate(x**2, x)} + C")
    print(f"  ∫ 1/x dx = {integrate(1/x, x)} + C")
    print(f"  ∫ e^x dx = {integrate(exp(x), x)} + C")
    print(f"  ∫ sin(x) dx = {integrate(sin(x), x)} + C")
    print(f"  ∫ cos(x) dx = {integrate(cos(x), x)} + C")
    
    # More complex integrals
    print(f"\nMore complex:")
    print(f"  ∫ x * e^x dx = {integrate(x * exp(x), x)} + C")
    print(f"  ∫ x^2 * sin(x) dx = {integrate(x**2 * sin(x), x)} + C")
    
    # Definite integrals
    print(f"\nDefinite integrals:")
    print(f"  ∫₀¹ x^2 dx = {integrate(x**2, (x, 0, 1))}")
    print(f"  ∫₀^π sin(x) dx = {integrate(sin(x), (x, 0, pi))}")
    print(f"  ∫₀^∞ e^(-x) dx = {integrate(exp(-x), (x, 0, oo))}")
    
    # Gaussian integral
    print(f"\nGaussian integral:")
    print(f"  ∫₋∞^∞ e^(-x^2) dx = {integrate(exp(-x**2), (x, -oo, oo))}")
    
    # Improper integral
    print(f"\nImproper integrals:")
    print(f"  ∫₁^∞ 1/x^2 dx = {integrate(1/x**2, (x, 1, oo))}")
    
    return integrate(x**2, x)


# ==============================================================================
# Example 5: Limits and Series
# ==============================================================================

def example_limits_series():
    """Limits and Taylor series expansions."""
    print("\n" + "=" * 60)
    print("Example 5: Limits and Series")
    print("=" * 60)
    
    x = symbols('x')
    
    # Limits
    print(f"\nLimits:")
    print(f"  lim(x→0) sin(x)/x = {limit(sin(x)/x, x, 0)}")
    print(f"  lim(x→0) (1-cos(x))/x^2 = {limit((1-cos(x))/x**2, x, 0)}")
    print(f"  lim(x→∞) (1 + 1/x)^x = {limit((1 + 1/x)**x, x, oo)}")
    print(f"  lim(x→0+) x*ln(x) = {limit(x*log(x), x, 0, '+')}")
    
    # One-sided limits
    print(f"\nOne-sided limits:")
    print(f"  lim(x→0+) 1/x = {limit(1/x, x, 0, '+')}")
    print(f"  lim(x→0-) 1/x = {limit(1/x, x, 0, '-')}")
    
    # Derivative definition via limit
    print(f"\nDerivative as limit:")
    h = symbols('h')
    f = x**2
    deriv_limit = limit(((f.subs(x, x+h) - f) / h), h, 0)
    print(f"  f(x) = x^2")
    print(f"  f'(x) = lim(h→0) [f(x+h)-f(x)]/h = {deriv_limit}")
    
    # Taylor series
    print(f"\nTaylor series expansions:")
    print(f"  e^x ≈ {sp.series(exp(x), x, 0, 6)}")
    print(f"  sin(x) ≈ {sp.series(sin(x), x, 0, 8)}")
    print(f"  cos(x) ≈ {sp.series(cos(x), x, 0, 8)}")
    print(f"  ln(1+x) ≈ {sp.series(log(1+x), x, 0, 6)}")
    
    # Series at different point
    print(f"\nSeries at x=1:")
    print(f"  e^x ≈ {sp.series(exp(x), x, 1, 4)}")
    
    return limit(sin(x)/x, x, 0)


# ==============================================================================
# Example 6: Matrices and Linear Algebra
# ==============================================================================

def example_matrices():
    """Symbolic matrix operations."""
    print("\n" + "=" * 60)
    print("Example 6: Matrices and Linear Algebra")
    print("=" * 60)
    
    # Define symbolic matrix
    A = Matrix([[1, 2], [3, 4]])
    B = Matrix([[5, 6], [7, 8]])
    
    print(f"\nMatrix A:")
    print(A)
    print(f"\nMatrix B:")
    print(B)
    
    # Basic operations
    print(f"\nMatrix operations:")
    print(f"  A + B = {A + B}")
    print(f"  A * B = {A * B}")
    print(f"  3 * A = {3 * A}")
    print(f"  A^T = {A.T}")
    
    # Determinant and inverse
    print(f"\nDeterminant and inverse:")
    print(f"  det(A) = {A.det()}")
    print(f"  A^(-1) = {A.inv()}")
    print(f"  Verify: A * A^(-1) = {simplify(A * A.inv())}")
    
    # Eigenvalues and eigenvectors
    print(f"\nEigenvalues and eigenvectors:")
    eigenvals = A.eigenvals()
    print(f"  Eigenvalues: {eigenvals}")
    
    eigenvects = A.eigenvects()
    print(f"  Eigenvectors: {eigenvects}")
    
    # Solve linear system Ax = b
    b = Matrix([5, 6])
    x = A.solve(b)
    print(f"\nSolve Ax = b where b = {b.T}:")
    print(f"  x = {x}")
    print(f"  Verify: A*x = {A * x}")
    
    # Matrix with symbols
    a, b, c, d = symbols('a b c d')
    M = Matrix([[a, b], [c, d]])
    print(f"\nSymbolic matrix M = [[a,b],[c,d]]:")
    print(f"  det(M) = {M.det()}")
    print(f"  trace(M) = {M.trace()}")
    
    return A


# ==============================================================================
# Example 7: Differential Equations
# ==============================================================================

def example_differential_equations():
    """Solving ordinary differential equations."""
    print("\n" + "=" * 60)
    print("Example 7: Differential Equations")
    print("=" * 60)
    
    x = symbols('x')
    f = sp.Function('f')
    
    # First-order ODE: f'(x) = f(x)
    print(f"\nFirst-order ODE: f'(x) = f(x)")
    de1 = sp.Eq(f(x).diff(x), f(x))
    sol1 = sp.dsolve(de1, f(x))
    print(f"  General solution: {sol1}")
    
    # First-order ODE with initial condition
    print(f"\nWith initial condition f(0) = 1:")
    sol1_ic = sp.dsolve(de1, f(x), ics={f(0): 1})
    print(f"  Particular solution: {sol1_ic}")
    
    # Second-order ODE: f''(x) + f(x) = 0 (harmonic oscillator)
    print(f"\nSecond-order ODE: f''(x) + f(x) = 0")
    de2 = sp.Eq(f(x).diff(x, 2) + f(x), 0)
    sol2 = sp.dsolve(de2, f(x))
    print(f"  General solution: {sol2}")
    
    # Second-order ODE with initial conditions
    print(f"\nWith f(0) = 0, f'(0) = 1:")
    sol2_ic = sp.dsolve(de2, f(x), ics={f(0): 0, f(x).diff(x).subs(x, 0): 1})
    print(f"  Particular solution: {sol2_ic}")
    
    # Damped harmonic oscillator
    print(f"\nDamped harmonic oscillator: f''(x) + 2*f'(x) + 5*f(x) = 0")
    de3 = sp.Eq(f(x).diff(x, 2) + 2*f(x).diff(x) + 5*f(x), 0)
    sol3 = sp.dsolve(de3, f(x))
    print(f"  Solution: {sol3}")
    
    # Non-homogeneous ODE
    print(f"\nNon-homogeneous: f'(x) + 2*f(x) = sin(x)")
    de4 = sp.Eq(f(x).diff(x) + 2*f(x), sin(x))
    sol4 = sp.dsolve(de4, f(x))
    print(f"  Solution: {sol4}")
    
    return sol1


# ==============================================================================
# Example 8: Code Generation (Lambdify)
# ==============================================================================

def example_code_generation():
    """Convert symbolic expressions to executable functions."""
    print("\n" + "=" * 60)
    print("Example 8: Code Generation with Lambdify")
    print("=" * 60)
    
    x = symbols('x')
    
    # Define symbolic expression
    expr = x**2 * sin(x) + exp(-x**2/2)
    print(f"\nSymbolic expression:")
    print(f"  f(x) = x^2 * sin(x) + exp(-x^2/2)")
    
    # Convert to NumPy function
    f_numpy = sp.lambdify(x, expr, 'numpy')
    
    # Evaluate numerically
    x_vals = np.linspace(0, 10, 1000)
    y_vals = f_numpy(x_vals)
    
    print(f"\nNumerical evaluation:")
    print(f"  f(0) = {f_numpy(0):.6f}")
    print(f"  f(π) = {f_numpy(np.pi):.6f}")
    print(f"  f(1) = {f_numpy(1):.6f}")
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, y_vals, 'b-', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Lambdified Function: $f(x) = x^2 \\sin(x) + e^{-x^2/2}$')
    plt.grid(True, alpha=0.3)
    plt.savefig('sympy_lambdify_plot.png', dpi=150)
    print("\nSaved: sympy_lambdify_plot.png")
    plt.close()
    
    # Multiple variables
    y = symbols('y')
    expr2 = x**2 + y**3
    f2 = sp.lambdify((x, y), expr2, 'numpy')
    
    print(f"\nMultivariable function:")
    print(f"  g(x, y) = x^2 + y^3")
    print(f"  g(2, 3) = {f2(2, 3)}")
    
    # Generate C code
    from sympy.utilities.codegen import codegen
    [(c_name, c_code), (h_name, h_header)] = codegen(
        ('my_function', expr), 'C'
    )
    
    print(f"\nGenerated C code snippet:")
    print("-" * 40)
    print('\n'.join(c_code.split('\n')[:15]))
    print("-" * 40)
    
    # LaTeX output
    latex_str = sp.latex(expr)
    print(f"\nLaTeX representation:")
    print(f"  {latex_str}")
    
    return f_numpy


# ==============================================================================
# Main Execution
# ==============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("SymPy: Symbolic Mathematics Examples")
    print("=" * 70)
    
    try:
        example_expressions()
        example_solving()
        example_derivatives()
        example_integrals()
        example_limits_series()
        example_matrices()
        example_differential_equations()
        example_code_generation()
        
        print("\n" + "=" * 70)
        print("All examples completed successfully!")
        print("Generated plots:")
        print("  - sympy_lambdify_plot.png")
        print("=" * 70)
        
    except ImportError as e:
        print(f"\nImport Error: {e}")
        print("Please install required packages:")
        print("  pip install sympy numpy matplotlib")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
