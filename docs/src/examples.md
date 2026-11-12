## Examples

Examples of the setup procedure can be found in the scripts in the `Chides.jl/examples` directory.

### Book Chapter 2

The reproducible interpolation examples in
[`examples/chapter02`](https://github.com/appelo/Chides.jl/tree/main/examples/chapter02)
cover cardinal functions, coefficient scaling and translation, nonsmooth data,
smooth convergence, seminorm projection, and tensor interpolation. From the
checkout, see `examples/chapter02/README.md` for setup, numerical checks, headless
figure generation, and explicit export to the companion book. The examples have
a dedicated Julia environment and use the local Chides interpolation routines.

### Book Chapter 4

[`examples/chapter04`](https://github.com/appelo/Chides.jl/tree/main/examples/chapter04)
recreates the one-dimensional dissipative wave experiments using the reusable
`DissipativeWave1D` solver. It includes periodic sine-wave checks, Gaussian
convergence at CFL 0.8 and 1, numerical tests, and an optional comparison with
the original Fortran source. See its README for the unequal displacement and
velocity orders, boundary corrections, reference solution, and explicit book export.

`independent_modes.jl` reproduces the book's **Independent Mode Checks** using
Chides interpolation and the existing wave kernel. It advances all Fourier
coefficients, records advection and wave rates, and exports both figures and the
LaTeX rate definitions. Run `test_independent_modes.jl` for rational interpolation,
full-grid coefficient comparisons, and convergence checks. These examples require
Julia only; the chapter README gives the commands.

`step_comparison.jl` regenerates Figure 4.4: equidistant and Chebyshev Lagrange
interpolation and two-point Hermite interpolation of a step. It shares Chapter 2's
interpolation helpers, displays the full oscillation range, and uses a 256-bit
reference formula for the high-order Hermite curve. Its CSV diagnostics retain
the Float64 discrepancy separately. See `test_step_comparison.jl` and the chapter
README for verification and standalone export.

`square_pulse.jl` adds a companion to the archived long-time square-pulse
comparison. It evolves FD8/RK4 and Hermite orders 3, 23, 31 through 2,000 periods
using powers of the discrete Fourier updates. It documents the time-step choices
and corrected FD stencil, checks high-order arithmetic, and exports three newly
named figures without replacing any archived plot.

`pulse_convergence.jl` adds the square-pulse L1 convergence experiment at CFL
100/101, using the archived advection template's inclusive initialization and
step counts for `m=1:5` and ten meshes through 51,200 cells. It shares the discrete advection
propagator, uses FFTW, integrates the reconstructed error with quadrature refinement,
and generates a figure and rate table while preserving the archived results.
See `test_pulse_convergence.jl` and the Chapter 4 README for validation and export.

`pulse_order.jl` accompanies the archived Figure 4.9 by varying `m=1:15` on
100- and 1,000-cell meshes. It reuses the pulse convergence experiment, measures
continuous L1 error over each cell, and compares it with discrete L1 error
sampled at the grid nodes. The continuous plot shows the mass-error bound from
inclusive initialization; this bound does not apply to the nodal norm.
An additional figure integrates the error of periodic piecewise-linear
interpolation of the same nodal solution values. This uses exact cell integrals
against the square pulse and retains the continuous mass-error bound.
Higher-order matrix powers use 256-bit
arithmetic. `test_pulse_order.jl` checks maps, precision, and direct transfers;
export preserves the archived order plot.

`weno_comparison.jl` adds two continuous-corner transport examples to the book's
**Comparison to WENO** section. WENO5-JS/SSP(3,3) is compared with Hermite–Taylor
at `m=1,5`, using the archived grid sizes and explicit analytic initialization on
half-cell-shifted nodes. The figures show solutions, corner close-ups, dense
pointwise errors, and total variation. Hermite uses the existing interpolation
map, polynomial translation, and dense reconstruction; the new WENO code is
restricted to scalar positive-speed advection. The README specifies the unequal
coefficient counts and the WENO linear display reconstruction used for its
continuous error diagnostics. `test_weno_comparison.jl` checks conservation,
smooth convergence orders, all Hermite coefficients against Fourier evolution,
and refined diagnostics. Export preserves every archived WENO/Hermite figure.

### Book Chapter 8

`examples/chapter08` solves the periodic planar compressible Navier–Stokes
system from Section 8.2.1 using the reusable `CompressibleNavierStokes2D` solver.
It stores specific volume, two velocities, pressure, and their scaled spatial
derivatives. Both viscosity and heat conduction are included. Local RK4
substeps advance the tensor Hermite interpolants; interpolation and polynomial
products reuse Chides kernels. The example README explains the ideal-gas
closure, constant transport coefficients, acoustic/diffusive step estimates,
and sampled positivity checks.

Run `examples/chapter08/setup.jl` in its own project environment, then
`examples/chapter08/run.jl` for the decaying vortex. Use `--case manufactured`
for the forced viscous/thermal verification or `--case entropy` for exact Euler
transport. Use `--case double_shear` for two perturbed periodic shear layers,
or `examples/chapter08/movies.jl --case double_shear --frames 161` to generate
their vorticity/density movie on a refined display mesh.
`convergence.jl` refines meshes and local RK substeps and measures
continuous reconstruction errors. Outputs include dense field plots, integral
histories, CSV errors, and source metadata. The book's existing figures are not
modified. See `test/test_cns2d.jl` and `examples/chapter08/runtests.jl` for checks.
