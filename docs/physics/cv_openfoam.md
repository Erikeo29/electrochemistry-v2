# Modélisation CV : OpenFOAM

L'approche OpenFOAM repose sur la méthode des volumes finis (FVM), standard dans le domaine de la mécanique des fluides numérique (CFD).

## Solveur : `electroChemFoam`
Le solveur est une adaptation de `laplacianFoam` ou `scalarTransportFoam` pour inclure les conditions aux limites électrochimiques spécifiques.

## Particularités
- **Structure du cas** : Utilisation des répertoires standards `0/`, `constant/` et `system/`.
- **Propriétés de transport** : Le coefficient de diffusion $D$ est défini dans `constant/transportProperties`.
- **Conditions aux Limites (BC)** : Implémentation de la loi de Butler-Volmer via des `codedFixedValue` ou des bibliothèques externes (`echemFvPatchFields`).

## Avantages
OpenFOAM permet de passer facilement d'une simulation 2D à une simulation 3D complète et de coupler l'électrochimie avec l'hydrodynamique (électrochimie sous flux), ce qui est une extension naturelle de ce projet.
刻