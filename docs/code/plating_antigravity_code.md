# Code Source : Electroplating (Antigravity)

Extraits du solveur galvanostatique optimisé.

## 1. Équation de Laplace & Butler-Volmer
```python
# Potentiel dans l'électrolyte
F_form = sigma * inner(grad(phi_solve), grad(v)) * dx

# Surtension d'activation
eta = phi_solve - Constant(E_eq)

# Densité de courant de Butler-Volmer
j_bv = j0 * (exp(alpha_a * f * eta) - exp(-alpha_c * f * eta))

# Ajout du terme source à la surface de la cathode
F_form += j_bv * v * ds(cathode_marker)
```

## 2. Boucle de Bisection (Mode Galvanostatique)
```python
for iteration in range(max_iterations):
    V_anode_mid = (V_anode_min + V_anode_max) / 2
    V_anode_const.assign(V_anode_mid)
    
    solver.solve()
    j_obtained = compute_cathode_current()
    
    if j_obtained < j_target:
        V_anode_max = V_anode_mid
    else:
        V_anode_min = V_anode_mid
```

## 3. Calcul de l'Épaisseur (Loi de Faraday)
```python
# Facteur de conversion vers épaisseur (microns)
faraday_factor = (time_deposition * M_Ni * efficiency * 1e6) / (n * F_const * rho_Ni)
thickness = np.abs(j_cathode) * faraday_factor
```
