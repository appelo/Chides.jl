using Chides
using Plots
using TaylorSeries
using LinearAlgebra
pyplot()

include("chides_sim_1d_scalar.jl")

function setinitialdata!(u,x,m,h)
    # Use a Taylor seris + AD to compute the
    # initial data
    affine(a) = a + Taylor1(typeof(a),2*m+4)
    #

    #idfun(x) = sin(pi*x)
    idfun(x) = exp(-36*x^2)

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

function pdeupdate!(uhinterp::Array{Float64,1},uspacetime::Array{Float64,2},sim)
    m = sim.m
    q = sim.q
    uspacetime[1:2*m+2,1] .= uhinterp
    cflx  = sim.dt/sim.h
    for idt=1:q
        for idx=0:2*m+1-idt
            uspacetime[1+idx,1+idt] = cflx/idt*(idx+1)*uspacetime[1+idx+1,1+idt-1]
        end
    end
    # Add up the Taylor series for the function and m first derivatives
    # at the next half step.
    dfact = 1.0
    for l=1:q
        dfact *= 0.5
        for i = 1:m+1
            uhinterp[i] += dfact*uspacetime[i,1+l]
        end
    end
end

function main(nx,m)
#    nx = 11
#    m = 3
    XL = -1.0
    XR = 1.0
    TEND = 10.0
    CFL = 0.95

    t_sim = sim.Chides_simulation_1D(m,nx,XL,XR,TEND,CFL)
    setinitialdata!(t_sim.up,t_sim.xp,t_sim.m,t_sim.h)

    nref = 25
    uplot = zeros((nx-1)*nref+1)
    xplot = LinRange(XL,XR,(nx-1)*nref+1)

    uex = zeros((nx-1)*nref+1)
    primal2plot1d!(uex,t_sim.up,t_sim.Hmat,nref)

    anim = @animate for it = 1:t_sim.nt
        evolvep2d_1d_scalar!(t_sim.ud,t_sim.up,t_sim.utaylor,t_sim.uhinterp,t_sim.Hmat,t_sim.m,t_sim.nx,t_sim)
        evolved2p_1d_scalar!(t_sim.up,t_sim.ud,t_sim.utaylor,t_sim.uhinterp,t_sim.Hmat,t_sim.m,t_sim.nx,t_sim)
        if mod(it,1) == 0
             primal2plot1d!(uplot,t_sim.up,t_sim.Hmat,nref)
             plot(xplot,uplot,xlims=(XL,XR),ylims=(-0.1, 1.1))
        end
    end
    gif(anim, string(@__DIR__,"/transport_solution.gif"), fps = 6)
    primal2plot1d!(uplot,t_sim.up,t_sim.Hmat,nref)
    return maximum(abs.(uplot.-uex))
end

m = 5
nx = 6
println(main(nx,m))

# errors = zeros(6,5)
# for m = 1:6
#     for p = 1:5
#         nx = 2^p+1
#         err = main(nx,m)
#         #println([m, nx, err])
#         errors[m,p] = err
#     end
#     println(log2.(errors[m,1:end-1]./errors[m,2:end]))
#     println("\n")
# end
