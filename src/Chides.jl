module Chides

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


export hermitemap!
export taylor_to_hermite1d!

end
