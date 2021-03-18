#push!(LOAD_PATH,"../")
using Chides
using BenchmarkTools
using Plots
using TaylorSeries
pyplot()

include("apa.jl")

function setinitialdata!(u,x,m,h)
    # Use a Taylor seris + AD to compute the
    # initial data
    affine(a) = a + Taylor1(typeof(a),m+3)
    #
    idfun(x) = sin(pi*x)
    for j = 1:length(x)
        xidat = affine(x[j])
        Uidat = idfun(xidat)
        # When using Taylor1 the factorial is built in
        df = 1.0
        for i = 1:m+1
            u[i,j] = df*Uidat[i-1]
            df *= h
        end
    end
end

function evolvep2d!(ud::Array{Float64,2},up::Array{Float64,2},utaylor::Array{Float64,2},uhinterp::Array{Float64,1},
    Hmat::Array{Float64,2},m::Int64,nx::Int64)

    # Periodic BC
    utaylor[1:m+1,1] .= view(up,1:m+1,nx-1)
    utaylor[1:m+1,2] .= view(up,1:m+1,1)
    taylor_to_hermite1d!(view(uhinterp,1:2*m+2),view(utaylor,1:m+1,1:2),Hmat,m)
    ud[1:m+1,1] .= view(uhinterp,1:m+1)
    for i = 1:nx-1
        taylor_to_hermite1d!(view(uhinterp,1:2*m+2),view(up,1:m+1,i:i+1),Hmat,m)
        ud[1:m+1,i+1] .= view(uhinterp,1:m+1)
    end
    # # Periodic BC
    utaylor[1:m+1,1] .= view(up,1:m+1,nx)
    utaylor[1:m+1,2] .= view(up,1:m+1,2)
    taylor_to_hermite1d!(view(uhinterp,1:2*m+2),view(utaylor,1:m+1,1:2),Hmat,m)
    ud[1:m+1,nx+1] .= view(uhinterp,1:m+1)
    return nothing
end

function evolved2p!(up::Array{Float64,2},ud::Array{Float64,2},utaylor::Array{Float64,2},uhinterp::Array{Float64,1},
    Hmat::Array{Float64,2},m::Int64,nx::Int64)

    for i = 1:nx
        taylor_to_hermite1d!(view(uhinterp,1:2*m+2),view(ud,1:m+1,i:i+1),Hmat,m)
        #taylor_to_hermite1d!(uhinterp,up[1:m+1,i:i+1],Hmat,m)
        ud[1:m+1,i+1] .= view(uhinterp,1:m+1)
    end
    return nothing
end




function main()
    # Tests for one dimensional arrays
    # Array sizes

    nx = 101
    XL = -1.0
    XR = 1.0
    m = 5
    TEND = 10.0
    CFL = 0.8

    t_sim = apa.Chides_simulation_1D(m,nx,XL,XR,TEND,CFL)

    # # Set up a primal and a dual grid
    # x = LinRange(XL,XR,nx)
    # h = x[2]-x[1]
    #
    # # Hermite interpolation matrix
    # Hmat = zeros(Float64,2*(m+1),2*(m+1))
    # xl = -1.0
    # xr = 1.0
    # xc = 0.0
    # icase = 0
    # hermitemap!(Hmat,m,xl,xr,xc,icase)
    #
    # # Local grid functions
    # uhermite = zeros(Float64,2*(m+1))
    # utaylor = rand(Float64,m+1,2)
    #
    # # Global grid functions
    # # Think of grid indices x_i i = 0,...,nx
    # up = zeros(Float64,m+1,nx)
    # ud = zeros(Float64,m+1,nx+1)


    setinitialdata!(t_sim.up,t_sim.xp,t_sim.m,t_sim.h)



#    @code_warntype evolvep2d!(t_sim.ud,t_sim.up,t_sim.utaylor,t_sim.uhinterp,t_sim.Hmat,t_sim.m,t_sim.nx)

    # @btime evolvep2d!(t_sim.ud,t_sim.up,t_sim.utaylor,t_sim.uhinterp,t_sim.Hmat,t_sim.m,t_sim.nx)

    #ud = copy(t_sim.ud)
    #up = copy(t_sim.up)
    #utaylor = copy(t_sim.utaylor)
    #uhinterp = copy(t_sim.uhinterp)
    #Hmat = copy(t_sim.Hmat)

    # @code_warntype evolvep2d!(ud,up,utaylor,uhinterp,Hmat,m,nx)

    for i = 1:500
        #evolvep2d!(ud,up,utaylor,uhinterp,Hmat,m,nx)
        evolvep2d!(t_sim.ud,t_sim.up,t_sim.utaylor,t_sim.uhinterp,t_sim.Hmat,t_sim.m,t_sim.nx)
        evolvep2d!(t_sim.up,t_sim.ud,t_sim.utaylor,t_sim.uhinterp,t_sim.Hmat,t_sim.m,t_sim.nx)
    end
end

using Profile


@profile for i in 1:100 main() end
Juno.profiler()
Profile.clear()

#@code_warntype evolvep2d!(ud,up,utaylor,uhinterp,Hmat,m,nx)
