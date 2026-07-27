module Chides

using LinearAlgebra

const fortran_lib_str = Base.Filesystem.dirname(pathof(@__MODULE__)) * "/libchides.so"

@doc raw"""
    hermitemap!(hmap::Array{Float64,2}, m::Int64, xl::Float64, xr::Float64,
                xc::Float64, icase::Int64=0)

Compute the matrix that transfers scaled derivative data of order `0` through `m` at
`xl` and `xr` to scaled derivative data at `xc`.

With ``h = x_r - x_l``, the matrix ``T`` written into `hmap` satisfies, for every
polynomial ``p`` of degree ``2m+1``,

```math
h^k \frac{D^k p(x_c)}{k!} = \sum_{j=0}^{m} T_{kj}\, h^j \frac{D^j p(x_l)}{j!}
                          + T_{k,j+m+1}\, h^j \frac{D^j p(x_r)}{j!},
\qquad k = 0,\dots,2m+1 .
```

Because the data carries the factors ``h^j``, the matrix depends only on the position
of `xc` within `[xl, xr]` and not on the cell width. Solvers therefore build it once
on the reference cell `(-1, 1, 0)` and reuse it on every cell of the grid.

`hmap` is overwritten in place and must be preallocated with size `(2m+2, 2m+2)`.
This is a wrapper around the Fortran routine `hermite_` in `libchides.so`.

# Arguments
- `hmap::Array{Float64,2}`: `(2m+2, 2m+2)` output matrix, overwritten in place.
- `m::Int64`: number of derivatives carried at each of the two data points.
- `xl::Float64`: location of the left data.
- `xr::Float64`: location of the right data.
- `xc::Float64`: center of the interpolant. Ignored unless `icase == 0`.
- `icase::Int64=0`: collocate at `xl` (`icase < 0`), at an interior `xc` with
  `xl < xc < xr` (`icase == 0`), or at `xr` (`icase > 0`). In the two boundary cases
  the leading `(m+1, m+1)` block is the identity, since the data at that end is
  reproduced exactly.

# Examples
```jldoctest
julia> H = zeros(2, 2);

julia> hermitemap!(H, 0, -1.0, 1.0, 0.0, 0);

julia> H
2×2 Matrix{Float64}:
  0.5  0.5
 -1.0  1.0
```
"""
function hermitemap!(hmap::Array{Float64,2}, m::Int64, xl::Float64, xr::Float64, xc::Float64, icase::Int64=0)

    ccall((:hermite_, fortran_lib_str ), Cvoid, (Ref{Int64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Int64}), Ref(m),Ref(xl),Ref(xr),Ref(xc),hmap,Ref(icase))

end

"""
    taylor_to_hermite1d!(uhermite, utaylor, Hmat::Array{Float64,2}, m::Int64)

Apply the Hermite map `Hmat` to the scaled Taylor data of two adjacent grid points.

`utaylor` holds the `m+1` scaled derivatives at the left point in its first column and
at the right point in its second column. `uhermite` receives the `2m+2` coefficients of
the interpolating polynomial on the cell, that is `Hmat * [utaylor[:,1]; utaylor[:,2]]`,
and must have length `2m+2`. Both are overwritten in place.

The arguments are deliberately untyped so that `view`s into the larger grid arrays can
be passed without allocating. This is the inner loop of the solvers, so the product is
written as an explicit accumulation rather than a matrix multiply.

# Examples
```jldoctest
julia> H = zeros(4, 4); hermitemap!(H, 1, -1.0, 1.0, 0.0, 0);

julia> utaylor = [1.0 2.0; 0.5 0.5];   # value and scaled derivative at each end

julia> uhermite = zeros(4);

julia> taylor_to_hermite1d!(uhermite, utaylor, H, 1);

julia> uhermite ≈ H * [utaylor[:, 1]; utaylor[:, 2]]
true
```
"""
function taylor_to_hermite1d!(uhermite,utaylor,Hmat::Array{Float64,2},m::Int64)
  for j = 1:2*m+2
    @inbounds uhermite[j] = 0.0
  end
  for i = 1:m+1
    for j = 1:2*m+2
      @inbounds uhermite[j] += Hmat[j,i]*utaylor[i,1]+Hmat[j,m+1+i]*utaylor[i,2]
    end
  end
end

@doc raw"""
    primal2plot1d!(uplot::Array{Float64,1}, up::Array{Float64,2},
                   Hmat::Array{Float64,2}, nref::Int64)

Evaluate the piecewise Hermite interpolant of a 1D grid function on a refined grid.

`up` stores scaled Taylor data, `up[k+1, i]` ``= h^k D^k u(x_i) / k!`` for
``k = 0,\dots,m``, so that `m = size(up, 1) - 1` and `nx = size(up, 2)`. Each of the
`nx-1` cells is mapped to its degree-``2m+1`` interpolant with `Hmat` and sampled at
`nref` equispaced points, giving values at the cell endpoints and `nref-1` interior
points. `uplot` is overwritten in place and must have length `(nx-1) * nref + 1`.

`Hmat` must be the interior reference-cell map, as built by
`hermitemap!(Hmat, m, -1.0, 1.0, 0.0, 0)`.

The result is exact for polynomials of degree at most ``2m+1``; for smooth data the
pointwise error is ``O(h^{2m+2})``.

# Examples
```jldoctest
julia> m, nx, nref = 1, 3, 4;

julia> H = zeros(2m + 2, 2m + 2); hermitemap!(H, m, -1.0, 1.0, 0.0, 0);

julia> xp = range(-1, 1, length = nx); h = step(xp);

julia> up = zeros(m + 1, nx);   # data for u(x) = x^2

julia> for i = 1:nx
           up[1, i] = xp[i]^2
           up[2, i] = h * 2 * xp[i]
       end

julia> uplot = zeros((nx - 1) * nref + 1);

julia> primal2plot1d!(uplot, up, H, nref);

julia> xplot = range(-1, 1, length = (nx - 1) * nref + 1);

julia> maximum(abs.(uplot .- xplot .^ 2)) < 1e-14   # degree 2 <= 2m+1, so exact
true
```
"""
function primal2plot1d!(uplot::Array{Float64,1},up::Array{Float64,2},Hmat::Array{Float64,2},nref::Int64)

  m = size(up,1)-1
  ucof = zeros(Float64,2*m+2)
  cof2plot = zeros(nref-1,2*m+2)
  for j = 1:2*m+2
    for i = 1:nref-1
      cof2plot[i,j] = (-.5+1.0*i/nref)^(j-1)
    end
  end
  # Should be done in parallel...
  for i = 1:size(up,2)-1
    taylor_to_hermite1d!(ucof,view(up,1:m+1,i:i+1),Hmat,m)
    uplot[1+nref*(i-1)] = up[1,i]
    mul!(view(uplot,1+nref*(i-1)+1:1+nref*(i-1)+nref-1),cof2plot,ucof)
  end
  uplot[end] = up[1,end]
end

@doc raw"""
    primal2plot2d!(uplot::Array{Float64,2}, up::Array{Float64,4},
                   Hmat::Array{Float64,2}, nref::Int64)

Evaluate the piecewise Hermite interpolant of a 2D grid function on a refined grid.

`up` stores mixed scaled Taylor data, `up[i+1, j+1, ix, iy]`
``= h_x^i h_y^j \, \partial_x^i \partial_y^j u(x_{ix}, y_{iy}) / (i!\, j!)``, so that
`m = size(up, 1) - 1`, `nx = size(up, 3)` and `ny = size(up, 4)`. Each cell is
interpolated by two tensor-product passes of `Hmat`: first the ``y``-data at the left
and right of the cell is mapped to the midpoints of the bottom and top faces, then the
result is mapped to the cell center. `uplot` is overwritten in place and must have size
`((nx-1) * nref + 1, (ny-1) * nref + 1)`.

`Hmat` must be the interior reference-cell map, as built by
`hermitemap!(Hmat, m, -1.0, 1.0, 0.0, 0)`.

The result is exact for tensor-product polynomials of degree at most ``2m+1`` in each
variable; for smooth data the pointwise error is ``O(h^{2m+2})``.

# Examples
```jldoctest
julia> m, nx, nref = 1, 3, 4;

julia> H = zeros(2m + 2, 2m + 2); hermitemap!(H, m, -1.0, 1.0, 0.0, 0);

julia> xp = range(-1, 1, length = nx); h = step(xp);

julia> up = zeros(m + 1, m + 1, nx, nx);   # data for u(x,y) = x*y

julia> for j = 1:nx, i = 1:nx
           up[1, 1, i, j] = xp[i] * xp[j]
           up[2, 1, i, j] = h * xp[j]
           up[1, 2, i, j] = h * xp[i]
           up[2, 2, i, j] = h * h
       end

julia> np = (nx - 1) * nref + 1; uplot = zeros(np, np);

julia> primal2plot2d!(uplot, up, H, nref);

julia> xg = range(-1, 1, length = np);

julia> maximum(abs.(uplot .- [x * y for x in xg, y in xg])) < 1e-14
true
```
"""
function primal2plot2d!(uplot::Array{Float64,2},up::Array{Float64,4},Hmat::Array{Float64,2},nref::Int64)

   m = size(up,1)-1
   nx = size(up,3)
   ny = size(up,4)
   udata = zeros(2*m+2)
   t = zeros(2*m+2,2*m+2)
   smcofs = zeros(2*m+2,2*m+2)
   xloc = zeros(nref+1,nref+1)
   yloc = zeros(nref+1,nref+1)
   for j = 0:nref
      for i = 0:nref
         xloc[1+i,1+j] = (-0.5+i/nref)
         yloc[1+i,1+j] = (-0.5+j/nref)
      end
   end
   for j = 1:ny-1
      for i = 1:nx-1
         U = view(up,1:m+1,1:m+1,i:i+1,j:j+1)
         for idy = 0:m
            # First interpolate all y-data at the bottom and top of a cell
            # to the middle of the bottom and top face.
            view(udata,1:m+1) .= view(U,1:m+1,1+idy,1,1)
            view(udata,m+2:2*m+2) .= view(U,1:m+1,1+idy,2,1)
            mul!(view(smcofs,1:2*m+2,1+idy),Hmat,udata)
            view(udata,1:m+1) .= view(U,1:m+1,1+idy,1,2)
            view(udata,m+2:2*m+2) .= view(U,1:m+1,1+idy,2,2)
            mul!(view(smcofs,1:2*m+2,m+2+idy),Hmat,udata)
         end
         # Then interpolate all data to the middle of the cell
         for idx = 0:2*m+1
            view(udata,1:m+1)     .= view(smcofs,1+idx,1:m+1)
            view(udata,m+2:2*m+2) .= view(smcofs,1+idx,m+2:2*m+2)
            mul!(view(t,1+idx,1:2*m+2),Hmat,udata)
         end
         # This should be improved ...
         uploc = zeros(nref+1,nref+1)
         for idy = 0:2*m+1
            for idx = 0:2*m+1
               uploc .+= t[1+idx,1+idy].*xloc.^idx.*yloc.^idy
            end
         end
         uplot[1+(i-1)*nref:i*nref+1,1+(j-1)*nref:j*nref+1] .= uploc
      end
   end
end


export hermitemap!
export taylor_to_hermite1d!
export primal2plot1d!
export primal2plot2d!

end
