using Plots
using Chides
using TaylorSeries
using BenchmarkTools, Test,LinearAlgebra
pyplot()


function d2p_2d!(up,ud,m,nx,ny,nref,Hmat,sim)
   q = sim[5]
   udata = zeros(2*m+2)
   smcofs = zeros(2*m+2,2*m+2)
   hip = zeros(2*m+2,2*m+2)
   uspacetime = zeros(2*m+2,2*m+2,q+1)
    for j = 1:ny
      for i = 1:nx
         @inbounds UD = view(ud,1:m+1,1:m+1,i:i+1,j:j+1)
         t2h_2d!(hip,UD,m,Hmat,udata,smcofs)
         pdeupdate!(hip,uspacetime,sim)
         @inbounds view(up,1:m+1,1:m+1,i,j) .= view(hip,1:m+1,1:m+1)
      end
   end
end

function p2d_2d!(ud,up,m,nx,ny,nref,Hmat,sim)
   q = sim[5]
   udata = zeros(2*m+2)
   smcofs = zeros(2*m+2,2*m+2)
   hip = zeros(2*m+2,2*m+2)
   uspacetime = zeros(2*m+2,2*m+2,q+1)
   UC = zeros(m+1,m+1,2,2)
    for j = 1:ny-1
      for i = 1:nx-1
         UP = view(up,1:m+1,1:m+1,i:i+1,j:j+1)
         UD = view(ud,1:m+1,1:m+1,i+1,j+1)
         t2h_2d!(hip,UP,m,Hmat,udata,smcofs)
         pdeupdate!(hip,uspacetime,sim)
         UD .= view(hip,1:m+1,1:m+1)
      end
   end
   # Left side
   i = 1
   for j = 1:ny-1
      UC[1:m+1,1:m+1,1,1:2] .= view(up,1:m+1,1:m+1,nx-1,j:j+1)
      UC[1:m+1,1:m+1,2,1:2] .= view(up,1:m+1,1:m+1,i,j:j+1)
      UD = view(ud,1:m+1,1:m+1,1,j+1)
      t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
      pdeupdate!(hip,uspacetime,sim)
      UD .= view(hip,1:m+1,1:m+1)
   end
   # Right side
   i = nx
   for j = 1:ny-1
      UC[1:m+1,1:m+1,1,1:2] .= view(up,1:m+1,1:m+1,i,j:j+1)
      UC[1:m+1,1:m+1,2,1:2] .= view(up,1:m+1,1:m+1,2,j:j+1)
      UD = view(ud,1:m+1,1:m+1,nx+1,j+1)
      t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
      pdeupdate!(hip,uspacetime,sim)
      UD .= view(hip,1:m+1,1:m+1)
   end
   # Bottom side
   j = 1
   for i = 1:nx-1
      UC[1:m+1,1:m+1,1:2,1] .= view(up,1:m+1,1:m+1,i:i+1,ny-1)
      UC[1:m+1,1:m+1,1:2,2] .= view(up,1:m+1,1:m+1,i:i+1,j)
      UD = view(ud,1:m+1,1:m+1,i+1,1)
      t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
      pdeupdate!(hip,uspacetime,sim)
      UD .= view(hip,1:m+1,1:m+1)
   end
   # Top side
   j = ny
   for i = 1:nx-1
      UC[1:m+1,1:m+1,1:2,1] .= view(up,1:m+1,1:m+1,i:i+1,j)
      UC[1:m+1,1:m+1,1:2,2] .= view(up,1:m+1,1:m+1,i:i+1,2)
      #UD = view(ud,1:m+1,1:m+1,i+1,ny+1)
      t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
      pdeupdate!(hip,uspacetime,sim)
      view(ud,1:m+1,1:m+1,i+1,ny+1) .= view(hip,1:m+1,1:m+1)
   end
   # Corners low-left
   UC[1:m+1,1:m+1,1,1] .= view(up,1:m+1,1:m+1,nx-1,ny-1)
   UC[1:m+1,1:m+1,2,1] .= view(up,1:m+1,1:m+1,1,ny-1)
   UC[1:m+1,1:m+1,1,2] .= view(up,1:m+1,1:m+1,nx-1,1)
   UC[1:m+1,1:m+1,2,2] .= view(up,1:m+1,1:m+1,1,1)
   t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
   pdeupdate!(hip,uspacetime,sim)
   UD = view(ud,1:m+1,1:m+1,1,1)
   UD .= view(hip,1:m+1,1:m+1)

   # Corners low-right
   UC[1:m+1,1:m+1,1,1] .= view(up,1:m+1,1:m+1,nx,ny-1)
   UC[1:m+1,1:m+1,2,1] .= view(up,1:m+1,1:m+1,2,ny-1)
   UC[1:m+1,1:m+1,1,2] .= view(up,1:m+1,1:m+1,nx,1)
   UC[1:m+1,1:m+1,2,2] .= view(up,1:m+1,1:m+1,2,1)
   t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
   pdeupdate!(hip,uspacetime,sim)
   UD = view(ud,1:m+1,1:m+1,nx+1,1)
   UD .= view(hip,1:m+1,1:m+1)

   # Corners top-left
   UC[1:m+1,1:m+1,1,1] .= view(up,1:m+1,1:m+1,nx-1,ny)
   UC[1:m+1,1:m+1,2,1] .= view(up,1:m+1,1:m+1,1,ny)
   UC[1:m+1,1:m+1,1,2] .= view(up,1:m+1,1:m+1,nx-1,2)
   UC[1:m+1,1:m+1,2,2] .= view(up,1:m+1,1:m+1,1,2)
   t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
   pdeupdate!(hip,uspacetime,sim)
   UD = view(ud,1:m+1,1:m+1,1,ny+1)
   UD .= view(hip,1:m+1,1:m+1)

   # Corners top-right
   UC[1:m+1,1:m+1,1,1] .= view(up,1:m+1,1:m+1,nx,ny)
   UC[1:m+1,1:m+1,2,1] .= view(up,1:m+1,1:m+1,2,ny)
   UC[1:m+1,1:m+1,1,2] .= view(up,1:m+1,1:m+1,nx,2)
   UC[1:m+1,1:m+1,2,2] .= view(up,1:m+1,1:m+1,2,2)
   t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
   pdeupdate!(hip,uspacetime,sim)
   UD = view(ud,1:m+1,1:m+1,nx+1,ny+1)
   UD .= view(hip,1:m+1,1:m+1)

end

function t2h_2d!(hip,U,m,Hmat,udata,smcofs)

   @inbounds for idy = 0:m
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
   @inbounds for idx = 0:2*m+1
      view(udata,1:m+1)     .= view(smcofs,1+idx,1:m+1)
      view(udata,m+2:2*m+2) .= view(smcofs,1+idx,m+2:2*m+2)
      mul!(view(hip,1+idx,1:2*m+2),Hmat,udata)
   end
end


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
   for j = 1:ny
      for i = 1:nx
         x, y = set_variables("x y", order=3*m+3)
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

function get_uplot_2d!(uplot,up,m,nx,ny,nref,Hmat)
   for j = 1:ny-1
      for i = 1:nx-1
         U = view(up,1:m+1,1:m+1,i:i+1,j:j+1)
         udata = zeros(2*m+2)
         t = zeros(2*m+2,2*m+2)
         smcofs = zeros(2*m+2,2*m+2)
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
         # DEAA should be improved
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

nx = 11
ny = 11
m = 9

# u = [0.8876089187502644,0.37551366571505984,0.07359205690183224,0.008027980347946162,0.0010935558816456137,3.8310013117826e-5,1.925870345642622e-6,8.086506375626357e-8,3.288851080363031e-9]
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
get_uplot_2d!(uplotex,up,m,nx,ny,nref,Hmat)

anim = @animate for it = 1:nt
   p2d_2d!(ud,up,m,nx,ny,nref,Hmat,sim)
   d2p_2d!(up,ud,m,nx,ny,nref,Hmat,sim)

   get_uplot_2d!(uplot,up,m,nx,ny,nref,Hmat)
   plot(xplot,yplot,uplot,st=:surface,
   camera=(-30,30), xlims=(-1.0,1.0), ylims=(-1.0,1.0))
end

gif(anim, string(@__DIR__,"/transport_2d_solution.gif"), fps = 6)

plot(xplot,yplot,uplot.-uplotex,st=:surface,
camera=(-30,30), xlims=(-1.0,1.0), ylims=(-1.0,1.0))

println([hx,maximum(abs.(uplot.-uplotex))])
