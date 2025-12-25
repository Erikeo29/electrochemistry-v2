# Équations Clés : Électrodéposition (Plating)

La simulation de l'électrodéposition du Nickel (Ni) utilise un modèle de distribution de courant secondaire (Secondary Current Distribution).

## 1. Distribution du Potentiel
Dans l'électrolyte, la conservation de la charge suit la loi d'Ohm avec une conductivité $\sigma$ constante :

$$ \nabla \cdot (-\sigma \nabla \phi) = 0 $$

Où $\phi$ est le potentiel électrique dans l'électrolyte [V].

## 2. Cinétique à la Cathode (Butler-Volmer)
La densité de courant locale $j$ à la surface de la cathode est régie par l'équation de Butler-Volmer :

$$ j = j_0 \left[ \exp\left(\frac{\alpha_a n F}{RT} \eta\right) - \exp\left(-\frac{\alpha_c n F}{RT} \eta\right) \right] $$

Où :
- $j_0$ est la densité de courant d'échange [A/m²].
- $\eta = \phi - E_{eq}$ est la surtension d'activation.

## 3. Loi de Faraday (Épaisseur du dépôt)
L'épaisseur locale du dépôt $h$ après un temps $t$ est proportionnelle à la densité de courant locale et à l'efficacité faradique :

$$ h = \frac{j \cdot M \cdot t_{dep} \cdot \varepsilon}{n \cdot F \cdot \rho} $$

Où :
- $M$ est la masse molaire du Nickel (58.69 g/mol).
- $\rho$ est la masse volumique du Nickel (8908 kg/m³).
- $\varepsilon$ est l'efficacité cathodique ($\approx 0.98$).
- $t_{dep}$ est le temps de déposition [s].
