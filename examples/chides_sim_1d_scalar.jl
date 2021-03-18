module sim
using Chides
mutable struct Chides_simulation_1D
    m             ::Int64
    q             ::Int64
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
        q = 2*m+4
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

        new(m,q,nx,XL,XR,h,TEND,dt,CFL,nt,Hmat,up,ud,utaylor,uhinterp,xp,xd)
    end
end
end

function evolvep2d_1d_scalar!(ud::Array{Float64,2},up::Array{Float64,2},utaylor::Array{Float64,2},uhinterp::Array{Float64,1},
    Hmat::Array{Float64,2},m::Int64,nx::Int64,sim)

    q = sim.q
    uspacetime = zeros(2*m+2,q+1)
    # Periodic BC
    utaylor[1:m+1,1] .= view(up,1:m+1,nx-1)
    utaylor[1:m+1,2] .= view(up,1:m+1,1)
    taylor_to_hermite1d!(view(uhinterp,1:2*m+2),view(utaylor,1:m+1,1:2),Hmat,m)
    pdeupdate!(uhinterp,uspacetime,sim)
    ud[1:m+1,1] .= view(uhinterp,1:m+1)
    for i = 1:nx-1
        taylor_to_hermite1d!(view(uhinterp,1:2*m+2),view(up,1:m+1,i:i+1),Hmat,m)
        pdeupdate!(uhinterp,uspacetime,sim)
        ud[1:m+1,i+1] .= view(uhinterp,1:m+1)
    end
    # # Periodic BC
    utaylor[1:m+1,1] .= view(up,1:m+1,nx)
    utaylor[1:m+1,2] .= view(up,1:m+1,2)
    taylor_to_hermite1d!(view(uhinterp,1:2*m+2),view(utaylor,1:m+1,1:2),Hmat,m)
    pdeupdate!(uhinterp,uspacetime,sim)
    ud[1:m+1,nx+1] .= view(uhinterp,1:m+1)
    return nothing
end

function evolved2p_1d_scalar!(up::Array{Float64,2},ud::Array{Float64,2},utaylor::Array{Float64,2},uhinterp::Array{Float64,1},
    Hmat::Array{Float64,2},m::Int64,nx::Int64,sim)
    q = sim.q
    uspacetime = zeros(2*m+2,q+1)
    for i = 1:nx
        taylor_to_hermite1d!(view(uhinterp,1:2*m+2),view(ud,1:m+1,i:i+1),Hmat,m)
        pdeupdate!(uhinterp,uspacetime,sim)
        up[1:m+1,i] .= view(uhinterp,1:m+1)
    end
    return nothing
end
