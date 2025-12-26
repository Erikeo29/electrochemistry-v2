# Code Source : CV Firedrake (Python)

Ce code utilise la librairie **Firedrake** pour résoudre les équations de transport. Voici une structure plus complète permettant de reproduire la simulation.

## 1. Imports et Configuration
```python
import numpy as np
import os
import argparse
from firedrake import *  # Importation principale FEM

# Paramètres physiques par défaut
F = 96485.0      # Faraday [C/mol]
R = 8.314        # Gaz parfaits [J/mol/K]
T = 298.15       # Température [K]
```

## 2. Classe de Simulation
L'encapsulation dans une classe permet de gérer proprement l'état (temps, concentrations passées).

```python
class CVSimulation:
    def __init__(self, mesh_path):
        self.mesh = Mesh(mesh_path)
        
        # Espace fonctionnel P2 (Quadratique) pour plus de précision
        self.V = FunctionSpace(self.mesh, "CG", 2)
        
        # Fonctions Inconnues (Actuelles et Pas de temps précédent)
        self.c_R = Function(self.V, name="ferro")
        self.c_O = Function(self.V, name="ferri")
        self.c_R_n = Function(self.V)
        self.c_O_n = Function(self.V)
        
        # Fonctions Test
        self.v_R = TestFunction(self.V)
        self.v_O = TestFunction(self.V)

    def setup_solver(self):
        # Définition des conditions aux limites et du solveur non-linéaire
        pass
```

## 3. Implémentation de Butler-Volmer
Cette fonction calcule le flux à l'interface en fonction du potentiel appliqué $E(t)$.

```python
def butler_volmer_flux(self, E_t, c_R, c_O):
    # Surtension
    eta = E_t - self.E0
    f = self.n * F / (R * T)
    
    # Termes cinétiques
    k_ox = self.k0 * exp(self.alpha * f * eta)
    k_red = self.k0 * exp(-(1 - self.alpha) * f * eta)
    
    # Flux résultant (mol/m²/s)
    # Note : c_R et c_O sont les champs scalaires Firedrake
    flux = (k_ox * c_R - k_red * c_O)
    return flux
```

## 4. Boucle Temporelle (Time Loop)
C'est ici que la forme variationnelle est résolue à chaque pas de temps.

```python
def run(self, dt, n_steps):
    # Forme variationnelle (Backward Euler)
    F_R = ((self.c_R - self.c_R_n) / dt * self.v_R * dx
           + self.D * inner(grad(self.c_R), grad(self.v_R)) * dx
           + self.butler_volmer_flux(...) * self.v_R * ds(self.WE_ID))
           
    # Configuration du solveur de Newton
    problem = NonlinearVariationalProblem(F_R == 0, self.c_R, bcs=[...])
    solver = NonlinearVariationalSolver(problem, solver_parameters={
        "snes_type": "newtonls",
        "ksp_type": "preonly",
        "pc_type": "lu"  # Solveur direct efficace pour la 2D
    })
    
    for step in range(n_steps):
        t = step * dt
        E_t = self.get_potential(t)
        solver.solve()
        # Mise à jour pour le pas suivant
        self.c_R_n.assign(self.c_R)
```