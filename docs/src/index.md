# Chides.jl

CHIDES — the Charles Hermite Interpolation and Differential Equation Solver — implements
Hermite–Taylor methods for hyperbolic partial differential equations.

Grid functions are stored as scaled Taylor coefficients rather than point values: a 1D
field is held as ``h^k D^k u(x_i)/k!`` for ``k = 0, \dots, m``. Adjacent cell data is
combined into a degree-``2m+1`` Hermite interpolant, which is then advanced half a time
step onto a staggered dual grid; two half steps return to the original grid.

## Installation

The package wraps a Fortran shared library that **must be built before the package can
be loaded**:

```
cd src/chides_fortran
make
```

This produces `libchides.so` and moves it to `src/`, where `Chides` looks for it. The
makefile supplies the macOS SDK path automatically, so no extra flags are needed there.

If the link step nevertheless fails with `ld: library not found for -lSystem`, the
compiler is looking for an SDK that no longer exists; building with
`SDKROOT=$(xcrun --show-sdk-path)` in the environment resolves it.

## Getting started

The central object is the Hermite interpolation matrix, built once on a reference cell
and reused on every cell of the grid:

```julia
using Chides

m = 3                                    # derivatives carried per grid point
Hmat = zeros(2m + 2, 2m + 2)
hermitemap!(Hmat, m, -1.0, 1.0, 0.0, 0)
```

See [`hermitemap!`](@ref) for the definition of the matrix, and the scripts in
`Chides.jl/examples` for complete transport solvers in one and two dimensions.

## Contents

```@contents
```
