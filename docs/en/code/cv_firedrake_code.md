# Source Code: CV Firedrake (Python)

This code uses the **Firedrake** library to solve transport equations. Here is a more complete structure allowing you to reproduce the simulation.

## 1. Imports and Setup
```python
import numpy as np
import os
import argparse
from firedrake import *  # Main FEM import

# Physical Parameters
F = 96485.0      # Faraday [C/mol]
R = 8.314        # Ideal Gas Constant [J/mol/K]
T = 298.15       # Temperature [K]
```

## 2. Simulation Class
Encapsulating logic in a class ensures clean state management (time, history).

```python
class CVSimulation:
    def __init__(self, mesh_path):
        self.mesh = Mesh(mesh_path)
        
        # P2 (Quadratic) Function Space for higher accuracy
        self.V = FunctionSpace(self.mesh, "CG", 2)
        
        # Unknown Functions (Current and Previous Time Step)
        self.c_R = Function(self.V, name="ferro")
        self.c_O = Function(self.V, name="ferri")
        self.c_R_n = Function(self.V)
        self.c_O_n = Function(self.V)
        
        # Test Functions
        self.v_R = TestFunction(self.V)
        self.v_O = TestFunction(self.V)

    def setup_solver(self):
        # Define BCs and Non-linear solver
        pass
```

## 3. Butler-Volmer Implementation
Calculates the interfacial flux based on applied potential $E(t)$.

```python
def butler_volmer_flux(self, E_t, c_R, c_O):
    # Overpotential
    eta = E_t - self.E0
    f = self.n * F / (R * T)
    
    # Kinetic terms
    k_ox = self.k0 * exp(self.alpha * f * eta)
    k_red = self.k0 * exp(-(1 - self.alpha) * f * eta)
    
    # Resulting Flux (mol/m²/s)
    # Note: c_R and c_O are Firedrake scalar fields
    flux = (k_ox * c_R - k_red * c_O)
    return flux
```

## 4. Time Loop
Where the variational form is solved at each time step.

```python
def run(self, dt, n_steps):
    # Variational Form (Backward Euler)
    F_R = ((self.c_R - self.c_R_n) / dt * self.v_R * dx
           + self.D * inner(grad(self.c_R), grad(self.v_R)) * dx
           + self.butler_volmer_flux(...) * self.v_R * ds(self.WE_ID))
           
    # Newton Solver Configuration
    problem = NonlinearVariationalProblem(F_R == 0, self.c_R, bcs=[...])
    solver = NonlinearVariationalSolver(problem, solver_parameters={
        "snes_type": "newtonls",
        "ksp_type": "preonly",
        "pc_type": "lu"  # Direct solver, efficient for 2D
    })
    
    for step in range(n_steps):
        t = step * dt
        E_t = self.get_potential(t)
        solver.solve()
        # Update for next step
        self.c_R_n.assign(self.c_R)
```
