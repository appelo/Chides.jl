# CHIDES - The Charles Hermite Interpolation and Differential Equation Solver 

Hermite–Taylor methods for hyperbolic partial differential equations. Grid functions are
stored as scaled Taylor coefficients rather than point values; adjacent cell data is
combined into a high-order Hermite interpolant and advanced half a time step onto a
staggered dual grid.

## Installation

Chides is not registered, so clone the repository and use it in place:

```bash
git clone https://github.com/appelo/Chides.jl.git
cd Chides.jl
```

The Julia code wraps a Fortran shared library, which **must be built before the package
can be used**:

```bash
cd src/chides_fortran
make
cd ../..
```

This compiles `libchides.so` and moves it to `src/`, where `Chides` looks for it. The
makefile supplies the macOS SDK path automatically, so no extra flags are needed there.
Rebuild after any change to the `.f90` sources; `make clean` removes the build products.

You can then load the package from the repository:

```bash
julia --project=.
```
```julia
julia> using Chides

julia> Hmat = zeros(2, 2);                        # m = 0: function values only

julia> hermitemap!(Hmat, 0, -1.0, 1.0, 0.0, 0);   # map data at -1 and 1 to 0

julia> Hmat
2×2 Matrix{Float64}:
  0.5  0.5
 -1.0  1.0
```

For general `m`, `Hmat` is `(2m+2)×(2m+2)` and carries `m` derivatives at each end;
see the documentation for `hermitemap!`.

## Dissipative wave examples

`DissipativeWave1D` solves `u_t=v`, `v_t=c^2*u_xx` in one dimension with
periodic or homogeneous Dirichlet boundaries. It reuses the interpolation
library and evolves the local polynomials in Julia:

```julia
s = DissipativeWave1D(3, 40; cfl=0.8)
initialize_wave!(s, x -> sin(2pi*x), x -> 2pi*cos(2pi*x))
solve_wave!(s; nsteps=50)
```

See [Chapter 4 examples](examples/chapter04/README.md) for the unequal displacement
and velocity orders, Gaussian convergence studies, Fortran comparisons, and
book figure generation. The same directory's `independent_modes.jl` reproduces
the book's periodic advection and wave mode checks, exporting figures and measured
rates from Julia. See its README for the standalone command and numerical tests.
`step_comparison.jl` in the same directory generates the three panels of Figure
4.4, including high-precision Hermite evaluation and interpolation diagnostics.

## Testing

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

If a `startup.jl` in your Julia configuration interferes with the test subprocess — for
example one that loads Revise, which can fail to precompile there and abort the run
before any test executes — disable it:

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test()'
```

The docstring examples in `src/Chides.jl` are executed as doctests during a
documentation build, so building the docs (below) exercises more of the library than the
test suite does.

## Documentation

The documentation is not yet published, so build it locally. This requires the Fortran
library to be built first, because the doctests call into it:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The first command only needs to be run once. `Pkg.develop` is required rather than
optional, since Chides is unregistered and cannot be resolved from a registry. The build
writes plain HTML to `docs/build/`, which opens directly in a browser:

```bash
open docs/build/index.html          # macOS; use xdg-open on Linux
```

Docstrings are also available in the REPL without building anything:

```julia
julia> using Chides

help?> hermitemap!
```

## Examples

Book examples in [`examples/chapter04`](examples/chapter04/README.md) include
wave convergence, independent mode checks, step interpolation, and square-pulse
comparisons of long-time behavior and convergence. The book's archived plots
and rates are preserved.

Complete transport solvers in one and two dimensions are in `examples/`. These additionally
require `Plots`, `TaylorSeries` and `BenchmarkTools`, which are not dependencies of the
package, and several of them write animated `.gif` output next to the script.

`old_rst_stuff/` holds the methodology write-ups for the original Fortran CHIDES code and
remains the best reference for the underlying mathematics.

### Periodic compressible Navier–Stokes

`CompressibleNavierStokes2D` implements the planar viscous and heat-conducting
system from book Section 8.2.1 in `(r,u,v,p)`, with `r=1/rho`. Initialize smooth
periodic data with `initialize_cns!`, then use `solve_cns!` or `step_cns!`.
The [Chapter 8 examples](examples/chapter08/README.md) provide a decaying vortex,
manufactured-solution and entropy-wave checks, a periodic double shear layer,
dense error measurements, and vorticity/density movies. Reproducible setup/run
commands accompany the examples. The method uses local Hermite–RK4 evolution
and constant transport coefficients; it does not include shock or wall models.
