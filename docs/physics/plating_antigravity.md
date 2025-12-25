# Modélisation Électrodéposition : Python (Antigravity)

Ce module simule la croissance d'une couche de Nickel sur un substrat complexe.

## Distribution Secondaire de Courant
Le modèle résout l'équation de Laplace pour le potentiel électrique. La cinétique de Butler-Volmer impose une relation non-linéaire entre le potentiel local et la densité de courant.

## Algorithme Galvanostatique (Recherche de DDC)
Pour correspondre aux conditions expérimentales où l'on impose un courant total (ex: 8 A/dm²), le code utilise une méthode de **bissection optimisée** :
1. Estimation d'un potentiel d'anode.
2. Résolution du champ de potentiel.
3. Intégration du courant sur la cathode.
4. Ajustement itératif du potentiel d'anode jusqu'à atteindre le courant cible.

## Résultats de Simulation
- **Carte d'épaisseur** : Visualisation 3D (exportable en VTK) de l'épaisseur déposée.
- **Uniformité** : Calcul de l'écart-type de l'épaisseur sur toute la surface.
- **Histogramme** : Analyse statistique de la répartition du métal.
