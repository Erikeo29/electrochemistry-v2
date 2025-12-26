# Code Source : Electroplating (Antigravity)

Ce code résout la **distribution de courant secondaire** pour prédire l'épaisseur de dépôt.

## 1. Imports et Maillage
Le chargement du maillage (souvent `.msh` de GMSH) est la première étape critique.

```python
from firedrake import *
import numpy as np

# Chargement du maillage avec les étiquettes physiques (markers)
mesh = Mesh("electrode_3d.msh")
V = FunctionSpace(mesh, "CG", 1) # Éléments linéaires P1 pour le potentiel

# Identifiants des surfaces (définis dans GMSH)
ID_CATHODE = 10
ID_ANODE = 11
```

## 2. Définition du Problème Variationnel
On cherche le potentiel $\phi$ tel que $\nabla \cdot (\sigma \nabla \phi) = 0$.

```python
# Inconnue et fonction test
phi = Function(V)
v = TestFunction(V)

# Paramètres
sigma = Constant(20.0) # Conductivité [S/m]
E_eq = Constant(-0.26) # Potentiel équilibre Nickel

# Forme bilinéaire (Loi d'Ohm dans le volume)
F = sigma * inner(grad(phi), grad(v)) * dx

# Ajout du terme non-linéaire de surface (Butler-Volmer)
# j(phi) est intégré sur la surface de la cathode
eta = phi - E_eq
j_bv = j0 * (exp(alpha_a * f * eta) - exp(-alpha_c * f * eta))
F += j_bv * v * ds(ID_CATHODE)

# Condition de Dirichlet à l'anode (Potentiel imposé)
bc_anode = DirichletBC(V, Constant(0.0), ID_ANODE)
```

## 3. Algorithme Galvanostatique (Recherche de $V_{anode}$)
On veut un courant total $I_{target}$ précis. Comme le problème est non-linéaire monotone, une méthode de bissection est robuste.

```python
def solve_galvanostatic(target_current_density):
    V_min, V_max = -2.0, 2.0
    
    for i in range(50):
        V_mid = (V_min + V_max) / 2
        
        # Mettre à jour la BC anode
        bc_anode.function_arg = Constant(V_mid)
        
        # Résoudre
        solve(F == 0, phi, bcs=[bc_anode])
        
        # Calculer le courant obtenu
        current_obtained = assemble(j_bv * ds(ID_CATHODE))
        
        # Ajuster les bornes
        if current_obtained < target_current_density:
            V_min = V_mid # Besoin de plus de potentiel
        else:
            V_max = V_mid
            
    return phi
```