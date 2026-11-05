# Chapter 2: Hermite interpolation

These scripts accompany Chapter 2 of *Hermite Interpolation and Differential
Equations*. They call the local Chides Fortran-backed interpolation routines;
independent formulas and exact polynomial calculations check their results.

## Setup and execution

Install Julia (the recorded environment uses Julia 1.13.3), GNU Fortran, and Make.
From the Chides checkout:

```sh
make -C src/chides_fortran
julia --startup-file=no --project=examples/chapter02 examples/chapter02/setup.jl
julia --startup-file=no --project=examples/chapter02 examples/chapter02/runtests.jl
julia --startup-file=no --project=examples/chapter02 examples/chapter02/run_all.jl
```

The dedicated `Project.toml` and resolved `Manifest.toml` isolate the example
packages. Setup registers this checkout with a relative path and instantiates its
dependencies. A different Julia version can require fresh resolution. The first
plotting run may take longer while Plots/GR compiles. Rendering uses GR's headless
mode and requires no display server.

Run any script separately with the same `--project` argument. Paths to helper
files and the default `output/` directory are independent of the working directory.
`--output /path/to/results` changes the output directory. Outputs contain five
PDF figures, a LaTeX tabular, CSV measurements, and `metadata.toml` recording
parameters, arithmetic, package versions, source commit, and working-tree status.
No timestamp is included; numeric outputs are deterministic on the recorded platform.

To export exactly the six manuscript assets, supply an explicit book root:

```sh
julia --startup-file=no --project=examples/chapter02 examples/chapter02/run_all.jl --book-root ..
```

The book root is resolved against the command's working directory. This command
replaces `hb_basis.pdf`, `hb_step.pdf`, `hb_abs.pdf`, `hb_convergence.pdf`,
`hb_tensor.pdf`, and `hb_projection_table.tex` in `Chapter_HB_definition/figs/`.
Other CSV and metadata files remain in the example output directory. Run the
numerical tests before exporting. Existing Python generators are historical
references and are not required by this workflow.

## Worked examples

| Script | Computation |
| --- | --- |
| `cardinal_basis.jl` | Degree-seven cardinal functions for unit scaled derivative data at each endpoint. |
| `scaling_and_translation.jl` | Polynomial reproduction on shifted cells of different widths; explicit triangular translation/rescaling, checked against independently computed Taylor coefficients. |
| `step_function.jl` | Two-point reconstruction of minus the sign function for m=1,5,15; compare degree-matched equispaced Lagrange interpolation. Show full excursions and a separate jump detail. |
| `absolute_value.jl` | Chides curves for m=1,2,4,8,16; compare cusp values with the exact central-binomial formula through m=512. |
| `smooth_convergence.jl` | Value and first-derivative convergence for sin(2 pi x), m=0,...,3 and 2,...,32 cells. |
| `projection_identity.jl` | Orthogonality, Pythagorean identity, and minimality checks for m=0,...,4, with an independent exact rational construction. |
| `tensor_interpolation.jl` | Reconstruct exp(x) cos(pi y) on an 8-by-6 rectangular cell grid, and measure the effect of omitting mixed derivatives. |

## Coefficients and independent checks

At a cell endpoint, store `h^k f^(k)/k!`. `Taylor1` supplies the factorial
normalization; multiply its coefficient of degree k by `h^k` only. The centered
Chides coordinate is `(x-center)/h`, whereas the chapter's endpoint coordinate is
`(x-a)/h`. For tensor interpolation use `hx^k hy^l/(k! l!)` and retain all mixed
partials through order m in each variable. The tensor example uses analytic
mixed derivatives; its independent reproduction test uses `x^5 y^4` on unequal
cell widths.

The numerical suite checks degrees through `2m+1`, center changes, all cubic
midpoint coefficients, endpoint matrix blocks, right-cardinal signs, sampled
total variation, exact reference curves, projection identities, and tensor data.
Smooth errors are sampled at 65 points per cell, including endpoints. Error
values below `1e-12` (function) or `1e-11` (derivative) are flagged as potentially
roundoff-limited; rates use the last pair above these thresholds. They must be
within 0.25 of orders `2m+2` and `2m+1`, respectively. These are sampled errors,
not certified continuous supremum norms.

For the projection experiment, `f=x^(2m+2)+x^(2m+3)+x`, `p` is the reconstructed
polynomial, `e=f-p`, `D` is derivative order m+1, and `w=(x^2-1)^(m+1)`.
The four table columns report:

- `abs(<De,Dp>)/sqrt(||Df||² ||Dp||²)`;
- `abs(||Df||²-||Dp||²-||De||²)/||Df||²`;
- `abs(||D(p+w)||²-||Dp||²-||Dw||²)/max(||Df||²,||Dw||²)`;
- maximum coefficient discrepancy against the exact rational interpolant,
  divided by the larger of one and its largest coefficient magnitude.

These integrals use 256-bit arithmetic on the Float64 reconstruction. The
independent rational calculation satisfies the three identities exactly.

The large-m cusp curve is an exact formula evaluated at high precision, **not**
a Chides computation through degree 1025. Finite-precision coefficient maps can
lose accuracy at large m. The displayed Chides curves stop at m=16 and report
their discrepancies against 256-bit reference values. The step and Lagrange
reconstructions use different data functionals; their comparison does not
establish stability of a nonlinear evolution scheme.
