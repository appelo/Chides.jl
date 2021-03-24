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
