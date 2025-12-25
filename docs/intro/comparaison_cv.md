# Comparaison des Modèles de Voltamétrie

Pour la simulation de la CV, deux approches numériques majeures ont été implémentées et comparées.

| Caractéristique | Python (Firedrake) | OpenFOAM (C++) |
| :--- | :--- | :--- |
| **Méthode** | Éléments Finis (FEM) | Volumes Finis (FVM) |
| **Précision** | Très élevée (Éléments P2) | Standard |
| **Flexibilité Géométrique** | Excellente (Maillage non-structuré) | Très bonne |
| **Vitesse de calcul** | Rapide pour la 2D | Optimisé pour la 3D massive |
| **Facilité d'édition** | Très simple (Python/DSL) | Plus complexe (C++/Compilation) |
| **Utilisation Type** | Recherche, prototypage rapide | Production, géométries industrielles |

## Analyse
Bien que les deux méthodes convergent vers les mêmes résultats physiques pour un maillage suffisamment fin, **Firedrake** est privilégié pour les études de précision cinétique grâce à ses éléments d'ordre élevé. **OpenFOAM** se distingue dès que l'on souhaite intégrer des effets de turbulence ou des écoulements fluides complexes autour de l'électrode.
