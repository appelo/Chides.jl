using Plots
using Chides
using TaylorSeries
using BenchmarkTools, Test,LinearAlgebra
gr()

include("transport_2d_funs.jl")

function pdeupdate!(uhinterp::Array{Float64,2},uspacetime::Array{Float64,3},sim)
   m = sim[4]
   q = sim[5]
   uspacetime .= 0.0
   uspacetime[1:2*m+2,1:2*m+2,1] .= uhinterp
   cflx  = sim[1]/sim[2]
   cfly  = sim[1]/sim[3]
   for idt=1:q
      for idy=0:2*m+1
        for idx=0:2*m
             uspacetime[1+idx,1+idy,1+idt] = cflx/idt*(idx+1.0)*uspacetime[1+idx+1,1+idy,1+idt-1]
          end
       end
      for idx=0:2*m+1
         for idy=0:2*m
            uspacetime[1+idx,1+idy,1+idt] = uspacetime[1+idx,1+idy,1+idt] + cfly/idt*(idy+1.0)*uspacetime[1+idx,1+idy+1,1+idt-1]
         end
      end
   end
   # Add up the Taylor series for the function and m first derivatives
   # at the next half step.
   dfact = 1.0
   for l=1:q
      dfact *= 0.5
      for j = 1:m+1
         for i = 1:m+1
            uhinterp[i,j] += dfact*uspacetime[i,j,1+l]
         end
      end
   end
end


function set_initial_data_2d!(up,xp,yp,hx,hy,nx,ny,m)
   # Use a Taylor seris + AD to compute the
   # interpolation data
   space = JetSpace(order = 3*m+3, variables = ["x", "y"])
   for j = 1:ny
      for i = 1:nx
         x, y = variables(space)
         x = x + xp[i]
         y = y + yp[j]
         xysq = exp(-36*(x*x+y*y)) #sin(pi*x+1.2)*sin(pi*y+2.3)
         dfy = 1.0
         for idy = 0:m
            dfx = 1.0
            for idx = 0:m
               up[1+idx,1+idy,i,j] = dfx*dfy*getcoeff(xysq, (idx,idy))
               dfx *=hx
            end
            dfy *=hy
         end
      end
   end
end


nx = 11
ny = 11
m = 4

up = zeros(Float64,m+1,m+1,nx,ny)
ud = zeros(Float64,m+1,m+1,nx+1,ny+1)

# Hermite interpolation matrix
Hmat = zeros(Float64,2*(m+1),2*(m+1))
xl = -1.0
xr = 1.0
xc = 0.0
icase = 0

hermitemap!(Hmat,m,xl,xr,xc,icase)

# Grids
XL = -1.0
XR =  1.0
YB = -1.0
YT =  1.0

xp = LinRange(XL,XR,nx)
yp = LinRange(YB,YT,ny)

hx = xp[2]-xp[1]
hy = yp[2]-yp[1]

set_initial_data_2d!(up,xp,yp,hx,hy,nx,ny,m)

# For plotting
nref = 20
z = LinRange(-0.5,0.5,nref+1)
xplot = zeros(nref*(nx-1)+1,nref*(ny-1)+1)
yplot = zeros(nref*(nx-1)+1,nref*(ny-1)+1)
uplot = zeros(nref*(nx-1)+1,nref*(ny-1)+1)

xloc = zeros(nref+1,nref+1)
yloc = zeros(nref+1,nref+1)
for j = 1:nref+1
  for i = 1:nref+1
    xloc[i,j] = z[i]
    yloc[i,j] = z[j]
    end
end
for j = 1:ny-1
   for i = 1:nx-1
      xplot[1+(i-1)*nref:i*nref+1,1+(j-1)*nref:j*nref+1] .= xp[i].+hx*(0.5.+xloc)
      yplot[1+(i-1)*nref:i*nref+1,1+(j-1)*nref:j*nref+1] .= yp[j].+hy*(0.5.+yloc)
   end
end


tend = 2.0
dt = 0.8*min(hx,hy)
nt = ceil(tend/dt)
dt = tend/nt

q = 4*m+4
sim = (dt,hx,hy,m,q)

uplotex = similar(uplot)
primal2plot2d!(uplotex,up,Hmat,nref)

anim = @animate for it = 1:nt
   p2d_2d!(ud,up,m,nx,ny,nref,Hmat,sim)
   d2p_2d!(up,ud,m,nx,ny,nref,Hmat,sim)

   primal2plot2d!(uplot,up,Hmat,nref)
   plot(xplot,yplot,uplot,st=:surface,
   camera=(-30,30), xlims=(-1.0,1.0), ylims=(-1.0,1.0))
end

gif(anim, string(@__DIR__,"/transport_2d_solution.gif"), fps = 6)

plot(xplot,yplot,uplot.-uplotex,st=:surface,
camera=(-30,30), xlims=(-1.0,1.0), ylims=(-1.0,1.0))

println([hx,maximum(abs.(uplot.-uplotex))])

primal2plot2d!(uplot,up,Hmat,nref)
pl = plot(xplot,yplot,uplot.-uplotex,st=:surface,camera=(-30,30), xlims=(-1.0,1.0), ylims=(-1.0,1.0))
display(pl)
