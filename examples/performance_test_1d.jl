using Chides
using BenchmarkTools

# Tests for one dimensional arrays

# Array sizes

nx = 100
m = 10

# Hermite interpolation matrix
Hmat = zeros(Float64,2*(m+1),2*(m+1))
xl = -1.0
xr = 1.0
xc = 0.0
icase = 0
hermitemap!(Hmat,m,xl,xr,xc,icase)

uhermite = zeros(Float64,2*(m+1))
utaylor = rand(Float64,m+1,2)
#utaylor[1,1] = 1.0
#utaylor[1,2] = 1.0


# Think of grid indices x_i i = 0,...,nx
up = zeros(Float64,m+1,nx+1)
ud = zeros(Float64,m+1,nx)



udata = reshape(utaylor,2*m+2,1)
uhermite1 = similar(uhermite)
uhermite2 = similar(uhermite)
uhermite3 = similar(uhermite)
uhermite4 = similar(uhermite)

@btime taylor_to_hermite1d1!(uhermite1,utaylor,Hmat,m)
@btime taylor_to_hermite1d2!(uhermite2,utaylor,Hmat,m)
@btime taylor_to_hermite1d3!(uhermite3,utaylor,Hmat,m)
@btime taylor_to_hermite1d4!(uhermite4,utaylor,Hmat,m)

println(isapprox(uhermite1,uhermite2))
println(isapprox(uhermite1,uhermite3))
println(isapprox(uhermite1,uhermite4))
println(isapprox(uhermite2,uhermite3))
println(isapprox(uhermite2,uhermite4))
println(isapprox(uhermite3,uhermite4))


include("apa.jl")

println("V1")
@btime interpolate_to_dual1d1!(ud,up,Hmat)
println("V2")
@btime interpolate_to_dual1d2!(ud,up,Hmat)
println("V10")
@btime interpolate_to_dual1d10!(ud,up,Hmat)

println("V4")
@btime interpolate_to_dual1d4!(ud,up,Hmat)
println("V5")
@btime interpolate_to_dual1d5!(ud,up,Hmat)
#println("V6")
#@btime interpolate_to_dual1d6!(ud,up,Hmat)

println("V7")
@btime interpolate_to_dual1d7!(ud,up,Hmat)
println("V8")
@btime interpolate_to_dual1d8!(ud,up,Hmat)
#println("V9")
#@btime interpolate_to_dual1d9!(ud,up,Hmat)
