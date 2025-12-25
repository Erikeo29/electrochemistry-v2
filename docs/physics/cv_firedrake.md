# Modélisation CV : Python & Firedrake

L'implémentation Firedrake utilise la puissance de la méthode des éléments finis (FEM) pour résoudre les équations de diffusion-réaction.

## Configuration du Modèle
- **Maillage** : Maillage 2D non-structuré généré avec GMSH, avec un raffinement local à l'interface électrode/électrolyte (échelle du micromètre).
- **Éléments** : Utilisation d'éléments Lagrange de degré 2 (P2) pour une meilleure précision sur les gradients de concentration.
- **Schéma temporel** : Euler implicite (Backward Euler) pour une stabilité inconditionnelle, crucial lors des balayages rapides.

## Résolution Numérique
Le système non-linéaire issu de la condition de Butler-Volmer est résolu à chaque pas de temps via un solveur de Newton-Raphson avec préconditionnement LU :
- `snes_type`: `newtonls`
- `pc_type`: `lu`

## Résultats Typiques
La simulation permet d'obtenir le profil de concentration du Ferrocyanure au cours du temps. On observe la formation de la couche de diffusion qui s'étend dans la solution "bulk".
Le voltammogramme résultant montre les pics caractéristiques anodiques et cathodiques, dont l'écart (Delta Ep) est un indicateur de la réversibilité cinétique.
