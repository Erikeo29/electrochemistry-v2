# Source Code: Electroplating (Antigravity)

This code solves the **secondary current distribution** to predict deposition thickness.

## 1. Imports and Mesh
Loading the mesh (often `.msh` from GMSH) is the first critical step.

```python
from firedrake import *
import numpy as np

# Load mesh with physical markers
mesh = Mesh("electrode_3d.msh")
V = FunctionSpace(mesh, "CG", 1) # Linear P1 elements for potential

# Surface IDs (defined in GMSH)
ID_CATHODE = 10
ID_ANODE = 11
```

## 2. Variational Problem Definition
We seek potential $\phi$ such that $\nabla · (\sigma  \nabla  \phi) = 0$.

```python
# Unknown and Test Function
phi = Function(V)
v = TestFunction(V)

# Parameters
sigma = Constant(20.0) # Conductivity [S/m]
E_eq = Constant(-0.26) # Equilibrium Potential Nickel

# Bilinear Form (Ohm's Law in volume)
F = sigma * inner(grad(phi), grad(v)) * dx

# Add non-linear surface term (Butler-Volmer)
# j(phi) is integrated over the cathode surface
eta = phi - E_eq
j_bv = j0 * (exp(alpha_a * f * eta) - exp(-alpha_c * f * eta))
F += j_bv * v * ds(ID_CATHODE)

# Dirichlet Condition at Anode (Imposed Potential)
bc_anode = DirichletBC(V, Constant(0.0), ID_ANODE)
```

## 3. Galvanostatic Algorithm (Finding $V_{anode}$)
We aim for a precise total current $I_{target}$. Since the problem is monotonic, a bisection method is robust.

```python
def solve_galvanostatic(target_current_density):
    V_min, V_max = -2.0, 2.0
    
    for i in range(50):
        V_mid = (V_min + V_max) / 2
        
        # Update Anode BC
        bc_anode.function_arg = Constant(V_mid)
        
        # Solve
        solve(F == 0, phi, bcs=[bc_anode])
        
        # Compute obtained current
        current_obtained = assemble(j_bv * ds(ID_CATHODE))
        
        # Adjust bounds
        if current_obtained < target_current_density:
            V_min = V_mid # Need more potential
        else:
            V_max = V_mid
            
    return phi
```
