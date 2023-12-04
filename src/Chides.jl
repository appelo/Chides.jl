module Chides

using LinearAlgebra

const fortran_lib_str = Base.Filesystem.dirname(pathof(@__MODULE__)) * "/libchides.so"

"""
 hermitemap!(hmap::Array{Float64,2}, m::Int64, xl::Float64, xr::Float64, xc::Float64, icase::Int64=0)

This subroutine computes the coefficient matrix tmat which transfers derivative data of order through m at xl and xr to the function and derivative data at xc. That is, if p is a polynomial of degree 2m+1

       h^k D^k p (xc)/k! = sum_(j=0)^m tmat(k,j) h^j D^j p(xl)/j!
                                 + tmat(k,j+m+1) h^j D^j p(xr)/j!

- icase < 0 => xc=xl (left boundary case)
- icase = 0 => xl < xc < xr
- icase > 0 => xc=xr (right boundary case)

# Arguments
- `hmap::Array{Float64,2}`: Hermite interpolation matrix
- `m::Int64`: Number of derivatives
- `xl::Float64`: Loccation of left data
- `xr::Float64`: Loccation of right data
- `xc::Float64`: Center of interpolant
- `icase::Int64=0`: Interpolation case, left, interior or right
"""
function hermitemap!(hmap::Array{Float64,2}, m::Int64, xl::Float64, xr::Float64, xc::Float64, icase::Int64=0)

    ccall((:hermite_, fortran_lib_str ), Cvoid, (Ref{Int64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Int64}), Ref(m),Ref(xl),Ref(xr),Ref(xc),hmap,Ref(icase))

end

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

function extrapolate!(p::Vector{Float64}, z::Float64, q::Int64)

    ccall((:extrapolate_, fortran_lib_str ), Cvoid, (Ref{Float64},Ref{Float64},Ref{Int64}), p,Ref(z),Ref(q))

end


export hermitemap!
export extrapolate!
export taylor_to_hermite1d!
export primal2plot1d!
export primal2plot2d!

end
