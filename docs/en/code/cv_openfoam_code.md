# Detailed Configuration: CV OpenFOAM

OpenFOAM simulation relies on precise dictionary configuration rather than scripting. Here are the essential files to reproduce the case.

## 1. Physical Properties (`constant/transportProperties`)
Defines diffusion coefficients and constants.
```cpp
/*--------------------------------*- C++ -*----------------------------------*\
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

// Diffusion Coefficient [m²/s]
D               7.0e-9;

// Number of electrons
n               1;

// Temperature [K]
T               298.15;
```

## 2. Butler-Volmer Boundary Condition (`0/C`)
To implement complex kinetics without recompiling a solver, `codedFixedValue` is used. It compiles a C++ snippet on the fly.

```cpp
boundaryField
{
    WE
    {
        type            codedFixedValue;
        value           uniform 0; // Init value
        name            butlerVolmerBC; // Unique name

        code
        #{
            const fvPatch& boundaryPatch = patch();
            const vectorField& Cf = boundaryPatch.Cf(); // Face centers
            scalarField& field = *this;

            // Get Time
            scalar t = this->db().time().value();
            
            // Parameters
            scalar E0 = 0.36;
            scalar k0 = 1e-5;
            scalar alpha = 0.5;
            scalar F = 96485.0;
            scalar R = 8.314;
            scalar Temp = 298.15;
            scalar f = F/(R*Temp);

            // Potential Signal (Triangular)
            scalar E_app = ...; // E(t) logic here

            // Flux calculation loop
            forAll(boundaryPatch, facei)
            {
                scalar eta = E_app - E0;
                // ... Butler-Volmer Calc ...
            }
        #};
    }
}
```

## 3. Mesh and Schemes (`system/`)
- `blockMeshDict`: Geometry generation.
- `fvSchemes`: Use `Euler` for time and `Gauss linear` for gradients (2nd order accuracy).

```
