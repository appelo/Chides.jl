using Plots
using Chides
using TaylorSeries
using BenchmarkTools, Test,LinearAlgebra
gr()

function set_initial_data_2d!(up,xp,yp,hx,hy,nx,ny,m)
   # Use a Taylor seris + AD to compute the
   # interpolation data
   space = JetSpace(order = 3*m+3, variables = ["x", "y"])
   for j = 1:ny
      for i = 1:nx
         x, y = variables(space)
         x = x + xp[i]
         y = y + yp[j]
         xysq = sin(pi*x+1.2)*sin(pi*y+2.3)
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

nx = 21
ny = 21
m = 3

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
nref = 10
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

function d2p_2d!(up,ud,m,nx,ny,nref,Hmat)
   udata = zeros(2*m+2)
   smcofs = zeros(2*m+2,2*m+2)
   hip = zeros(2*m+2,2*m+2)
    for j = 1:ny
      for i = 1:nx
         @inbounds UD = view(ud,1:m+1,1:m+1,i:i+1,j:j+1)
         t2h_2d!(hip,UD,m,Hmat,udata,smcofs)
         @inbounds view(up,1:m+1,1:m+1,i,j) .= view(hip,1:m+1,1:m+1)
      end
   end
end

function p2d_2d!(ud,up,m,nx,ny,nref,Hmat)
   udata = zeros(2*m+2)
   smcofs = zeros(2*m+2,2*m+2)
   UC = zeros(m+1,m+1,2,2)
   hip = zeros(2*m+2,2*m+2)
    for j = 1:ny-1
      for i = 1:nx-1
         UP = view(up,1:m+1,1:m+1,i:i+1,j:j+1)
         UD = view(ud,1:m+1,1:m+1,i+1,j+1)
         t2h_2d!(hip,UP,m,Hmat,udata,smcofs)
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
      UD .= view(hip,1:m+1,1:m+1)
   end
   # Right side
   i = nx
   for j = 1:ny-1
      UC[1:m+1,1:m+1,1,1:2] .= view(up,1:m+1,1:m+1,i,j:j+1)
      UC[1:m+1,1:m+1,2,1:2] .= view(up,1:m+1,1:m+1,2,j:j+1)
      UD = view(ud,1:m+1,1:m+1,nx+1,j+1)
      t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
      UD .= view(hip,1:m+1,1:m+1)
   end
   # Bottom side
   j = 1
   for i = 1:nx-1
      UC[1:m+1,1:m+1,1:2,1] .= view(up,1:m+1,1:m+1,i:i+1,ny-1)
      UC[1:m+1,1:m+1,1:2,2] .= view(up,1:m+1,1:m+1,i:i+1,j)
      UD = view(ud,1:m+1,1:m+1,i+1,1)
      t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
      UD .= view(hip,1:m+1,1:m+1)
   end
   # Top side
   j = ny
   for i = 1:nx-1
      UC[1:m+1,1:m+1,1:2,1] .= view(up,1:m+1,1:m+1,i:i+1,j)
      UC[1:m+1,1:m+1,1:2,2] .= view(up,1:m+1,1:m+1,i:i+1,2)
      #UD = view(ud,1:m+1,1:m+1,i+1,ny+1)
      t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
      view(ud,1:m+1,1:m+1,i+1,ny+1) .= view(hip,1:m+1,1:m+1)
   end
   # Corners low-left
   UC[1:m+1,1:m+1,1,1] .= view(up,1:m+1,1:m+1,nx-1,ny-1)
   UC[1:m+1,1:m+1,2,1] .= view(up,1:m+1,1:m+1,1,ny-1)
   UC[1:m+1,1:m+1,1,2] .= view(up,1:m+1,1:m+1,nx-1,1)
   UC[1:m+1,1:m+1,2,2] .= view(up,1:m+1,1:m+1,1,1)
   t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
   UD = view(ud,1:m+1,1:m+1,1,1)
   UD .= view(hip,1:m+1,1:m+1)

   # Corners low-right
   UC[1:m+1,1:m+1,1,1] .= view(up,1:m+1,1:m+1,nx,ny-1)
   UC[1:m+1,1:m+1,2,1] .= view(up,1:m+1,1:m+1,2,ny-1)
   UC[1:m+1,1:m+1,1,2] .= view(up,1:m+1,1:m+1,nx,1)
   UC[1:m+1,1:m+1,2,2] .= view(up,1:m+1,1:m+1,2,1)
   t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
   UD = view(ud,1:m+1,1:m+1,nx+1,1)
   UD .= view(hip,1:m+1,1:m+1)

   # Corners top-left
   UC[1:m+1,1:m+1,1,1] .= view(up,1:m+1,1:m+1,nx-1,ny)
   UC[1:m+1,1:m+1,2,1] .= view(up,1:m+1,1:m+1,1,ny)
   UC[1:m+1,1:m+1,1,2] .= view(up,1:m+1,1:m+1,nx-1,2)
   UC[1:m+1,1:m+1,2,2] .= view(up,1:m+1,1:m+1,1,2)
   t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
   UD = view(ud,1:m+1,1:m+1,1,ny+1)
   UD .= view(hip,1:m+1,1:m+1)

   # Corners top-right
   UC[1:m+1,1:m+1,1,1] .= view(up,1:m+1,1:m+1,nx,ny)
   UC[1:m+1,1:m+1,2,1] .= view(up,1:m+1,1:m+1,2,ny)
   UC[1:m+1,1:m+1,1,2] .= view(up,1:m+1,1:m+1,nx,2)
   UC[1:m+1,1:m+1,2,2] .= view(up,1:m+1,1:m+1,2,2)
   t2h_2d!(hip,UC,m,Hmat,udata,smcofs)
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

#@btime
p2d_2d!(ud,up,m,nx,ny,nref,Hmat)
@btime d2p_2d!(up,ud,m,nx,ny,nref,Hmat)

get_uplot_2d!(uplot,up,m,nx,ny,nref,Hmat)
#plot(xplot,yplot,uplot,st=:surface,camera=(-30,30))

plot(xplot,yplot,uplot.-sin.(pi*xplot.+1.2).*sin.(pi*yplot.+2.3),st=:surface,camera=(-30,30))
#
