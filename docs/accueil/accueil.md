# Bienvenue sur la Plateforme de Simulation Électrochimique

Cette application interactive regroupe des outils de modélisation avancés pour l'étude de phénomènes électrochimiques complexes. Elle a été conçue comme une référence académique et technique pour la simulation numérique.

## Objectifs de la Plateforme
- **Visualiser** les gradients de concentration et les distributions de potentiel.
- **Comparer** différentes approches numériques (Éléments Finis via Python/Firedrake vs Volumes Finis via OpenFOAM).
- **Analyser** l'impact des paramètres physiques sur les résultats expérimentaux (Voltammogrammes, Épaisseur de dépôt).

---

## 1. Voltamétrie Cyclique (CV)
La CV est une technique électroanalytique puissante. Nous modélisons ici le transport de masse instationnaire couplé à une réaction redox réversible.
- **Python (Firedrake)** : Utilisation de la méthode des éléments finis (FEM) pour une précision spatiale élevée.
- **OpenFOAM** : Approche orientée mécanique des fluides pour des géométries complexes.

## 2. Électrodéposition (Plating)
Simulation du dépôt électrolytique de Nickel sur des électrodes à géométrie spécifique (3-pad).
- **Mode Galvanostatique** : Recherche automatique du potentiel d'anode pour atteindre une densité de courant cible (DDC).
- **Uniformité** : Analyse de la distribution d'épaisseur et optimisation du procédé.

---
*Développé dans le cadre de projets de R&D - Version 1.0.0*
