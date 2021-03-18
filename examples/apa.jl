module apa
using Chides
mutable struct Chides_simulation_1D
    m            ::Int64
    nx            ::Int64
    XL            ::Float64
    XR            ::Float64
    h             ::Float64
    TEND          ::Float64
    dt            ::Float64
    CFL           ::Float64
    nt            ::Int64
    Hmat          ::Array{Float64,2}
    up            ::Array{Float64,2}
    ud            ::Array{Float64,2}
    utaylor       ::Array{Float64,2}
    uhinterp      ::Array{Float64,1}
    xp            ::Array{Float64,1}
    xd            ::Array{Float64,1}
    function Chides_simulation_1D(m::Int64,nx::Int64,
        XL::Float64,XR::Float64,TEND::Float64,CFL::Float64)


        Hmat = zeros(Float64,2*(m+1),2*(m+1))
        hermitemap!(Hmat,m,-1.0,1.0,0.0,0)

        xp = LinRange(XL,XR,nx)
        h = xp[2]-xp[1]
        xd = LinRange(XL-0.5*h,XR+0.5*h,nx+1)
        dt = CFL*h
        nt = ceil(TEND/dt)
        dt = TEND/nt
        uhinterp = zeros(Float64,2*(m+1))
        utaylor = zeros(Float64,m+1,2)
        # Global grid functions
        # Think of grid indices x_i i = 0,...,nx
        up = zeros(Float64,m+1,nx)
        ud = zeros(Float64,m+1,nx+1)

        new(m,nx,XL,XR,h,TEND,dt,CFL,nt,Hmat,up,ud,utaylor,uhinterp,xp,xd)
    end
end
end
