# Équations Clés : Voltamétrie Cyclique (CV)

La modélisation de la voltamétrie cyclique repose sur la résolution du transport de masse des espèces électroactives couplé à la cinétique de transfert de charge à l'interface électrode/électrolyte.

## 1. Transport de Masse
Dans une solution support en excès, la migration est négligée. Le transport est régi par la loi de diffusion de Fick :

$$ \frac{\partial c_i}{\partial t} = D_i \nabla^2 c_i $$

Où :
- $c_i$ est la concentration de l'espèce $i$ (Oxydant $O$ ou Réducteur $R$) [mol/m³].
- $D_i$ est le coefficient de diffusion [m²/s].

## 2. Cinétique à l'Électrode (Butler-Volmer)
Le flux molaire à la surface de l'électrode de travail (WE) est défini par la relation de Butler-Volmer :

$$ J = k_0 \left[ c_R \exp\left(\frac{\alpha n F}{RT} \eta\right) - c_O \exp\left(-\frac{(1-\alpha) n F}{RT} \eta\right) \right] $$

Avec la surtension $\eta = E(t) - E^0$.

### Paramètres de simulation (valeurs utilisées dans le projet) :
- Constante de Faraday $F = 96485$ C/mol
- Constante des gaz parfaits $R = 8.314$ J/(mol·K)
- Température $T = 298.15$ K
- Coefficient de transfert $\alpha = 0.5$
- Vitesse de réaction standard $k_0 \approx 10^{-5}$ m/s

## 3. Signal de Potentiel
Le potentiel appliqué $E(t)$ varie linéairement avec le temps (balayage triangulaire) :

$$ E(t) = \begin{cases} E_{start} + v \cdot t & \text{pour le balayage aller} \\ E_{vertex} - v \cdot t & \text{pour le balayage retour} \end{cases} $$

Où $v$ est la vitesse de balayage (scan rate) en V/s.
