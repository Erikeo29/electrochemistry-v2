# Configuration Détaillée : CV OpenFOAM

OpenFOAM ne nécessite pas de "code" au sens script, mais une configuration précise de dictionnaires. Voici les fichiers essentiels pour reproduire le cas.

## 1. Propriétés Physiques (`constant/transportProperties`)
C'est ici que l'on définit les coefficients de diffusion.
```cpp
/*--------------------------------	- C++ -*	----------------------------------*
| =========                 |                                                 |
| \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      transportProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Coefficient de diffusion [m²/s]
D               7.0e-9;

// Nombre d'électrons échangés
n               1;

// Température [K]
T               298.15;
```

## 2. Condition Limite Butler-Volmer (`0/C`)
Pour implémenter une cinétique complexe sans recompiler un solveur, on utilise `codedFixedValue`. Cela compile un petit morceau de C++ à la volée.

```cpp
boundaryField
{
    WE
    {
        type            codedFixedValue;
        value           uniform 0; // Valeur init
        name            butlerVolmerBC; // Nom unique

        code
        #{
            const fvPatch& boundaryPatch = patch();
            const vectorField& Cf = boundaryPatch.Cf(); // Centres des faces
            scalarField& field = *this;

            // Récupérer le temps
            scalar t = this->db().time().value();
            
            // Paramètres (Hardcodés ou lus depuis dictionnaire)
            scalar E0 = 0.36;
            scalar k0 = 1e-5;
            scalar alpha = 0.5;
            scalar F = 96485.0;
            scalar R = 8.314;
            scalar Temp = 298.15;
            scalar f = F/(R*Temp);

            // Signal de potentiel (Triangulaire)
            scalar E_app = ...; // Logique E(t) ici

            // Calcul du flux pour chaque face de l'électrode
            forAll(boundaryPatch, facei)
            {
                scalar eta = E_app - E0;
                // ... Calcul Butler-Volmer ...
                // Note: En FVM, on impose souvent le gradient ou une valeur mixte
                // field[facei] = ...; 
            }
        #};
    }
}
```

## 3. Maillage et Schémas (`system/`)
- `blockMeshDict` : Pour générer la géométrie (si simple).
- `fvSchemes` : Choisir `Euler` pour le temps et `Gauss linear` pour les gradients (ordre 2).

```