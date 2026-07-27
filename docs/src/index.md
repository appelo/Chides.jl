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

This produces `libchides.so` and moves it to `src/`, where `Chides` looks for it.

On macOS with a Homebrew `gfortran` the link step may fail with
`ld: library not found for -lSystem`, because the makefile passes no SDK path. In that
case link manually after `make` has compiled the object files:

```
gfortran -shared -L"$(xcrun --show-sdk-path)/usr/lib" -o libchides.so \
    polymul.o polyfuns.o Hermite.o
mv libchides.so ../
```

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
