# Configuration : CV OpenFOAM

Extraits des dictionnaires de configuration pour OpenFOAM.

## 1. Propriétés de Transport (`constant/transportProperties`)
```cpp
// Diffusion du Ferrocyanure
D           7.0e-9;
```

## 2. Condition de Butler-Volmer (`0/c_R`)
L'implémentation repose souvent sur un `codedFixedValue` pour le flux :
```cpp
patchName
{
    type            codedFixedValue;
    value           uniform 1.0;
    name            butlerVolmer;
    
    code
    #{
        // Calcul du flux en fonction de E(t)
        // J = k0 * (cR * exp(a*f*eta) - cO * exp(-ac*f*eta))
    #};
}
```

## 3. Schéma Temporel (`system/controlDict`)
```cpp
application     electroChemFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         100;
deltaT          0.01;
writeControl    runTime;
writeInterval   1;
```
