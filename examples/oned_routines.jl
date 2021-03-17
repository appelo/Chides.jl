function interpolate_to_dual1d1!(ud::Array{Float64,2},up::Array{Float64,2},Hmat::Array{Float64,2})
  m = size(ud,1)-1
  ucof = zeros(Float64,2*m+2)
  for i = 1:size(ud,2)
    taylor_to_hermite1d1!(ucof,up[1:m+1,i:i+1],Hmat,m)
    ud[1:m+1,i] .= ucof[1:m+1]
  end
end

function interpolate_to_dual1d2!(ud::Array{Float64,2},up::Array{Float64,2},Hmat::Array{Float64,2})
  m = size(ud,1)-1
  ucof = zeros(Float64,2*m+2)
  for i = 1:size(ud,2)
    taylor_to_hermite1d1!(ucof,view(up,1:m+1,i:i+1),Hmat,m)
    ud[1:m+1,i] .= ucof[1:m+1]
  end
end

function interpolate_to_dual1d3!(ud::Array{Float64,2},up::Array{Float64,2},Hmat::Array{Float64,2})
  m = size(ud,1)-1
  ucof = zeros(Float64,2*m+2)
  for i = 1:size(ud,2)
    taylor_to_hermite1d2!(ucof,up[1:m+1,i:i+1],Hmat,m)
    ud[1:m+1,i] .= ucof[1:m+1]
  end
end

function interpolate_to_dual1d4!(ud::Array{Float64,2},up::Array{Float64,2},Hmat::Array{Float64,2})
  m = size(ud,1)-1
  ucof = zeros(Float64,2*m+2)
  for i = 1:size(ud,2)
    taylor_to_hermite1d2!(ucof,view(up,1:m+1,i:i+1),Hmat,m)
    ud[1:m+1,i] .= ucof[1:m+1]
  end
end


function interpolate_to_dual1d5!(ud::Array{Float64,2},up::Array{Float64,2},Hmat::Array{Float64,2})
  m = size(ud,1)-1
  ucof = zeros(Float64,2*m+2)
  for i = 1:size(ud,2)
    taylor_to_hermite1d3!(ucof,up[1:m+1,i:i+1],Hmat,m)
    ud[1:m+1,i] .= ucof[1:m+1]
  end
end

function interpolate_to_dual1d6!(ud::Array{Float64,2},up::Array{Float64,2},Hmat::Array{Float64,2})
  m = size(ud,1)-1
  ucof = zeros(Float64,2*m+2)
  for i = 1:size(ud,2)
    taylor_to_hermite1d3!(ucof,view(up,1:m+1,i:i+1),Hmat,m)
    ud[1:m+1,i] .= ucof[1:m+1]
  end
end

function interpolate_to_dual1d7!(ud::Array{Float64,2},up::Array{Float64,2},Hmat::Array{Float64,2})
  m = size(ud,1)-1
  ucof = zeros(Float64,2*m+2)
  for i = 1:size(ud,2)
    taylor_to_hermite1d4!(ucof,view(up,1:m+1,i:i+1),Hmat,m)
    ud[1:m+1,i] .= ucof[1:m+1]
  end
end


function taylor_to_hermite1d1!(uhermite::Array{Float64,1},utaylor::Array{Float64,2},Hmat::Array{Float64,2},m::Int64)
  udata = similar(uhermite)
  udata[1:m+1] .= utaylor[1:m+1,1]
  udata[m+2:2*m+2] .= utaylor[1:m+1,2]
  uhermite .= Hmat*udata
end

function taylor_to_hermite1d1!(uhermite::Array{Float64,1},utaylor,Hmat::Array{Float64,2},m::Int64)
  udata = similar(uhermite)
  udata[1:m+1] .= utaylor[1:m+1,1]
  udata[m+2:2*m+2] .= utaylor[1:m+1,2]
  uhermite .= Hmat*udata
end

function taylor_to_hermite1d1!(uhermite,utaylor,Hmat::Array{Float64,2},m::Int64)
  udata = similar(uhermite)
  udata[1:m+1] .= utaylor[1:m+1,1]
  udata[m+2:2*m+2] .= utaylor[1:m+1,2]
  uhermite .= Hmat*udata
end

function taylor_to_hermite1d2!(uhermite,utaylor,Hmat::Array{Float64,2},m::Int64)
  uhermite .= view(Hmat,1:2*m+2,1:m+1)*view(utaylor,1:m+1,1)
  uhermite .+= view(Hmat,1:2*m+2,m+2:2*m+2)*view(utaylor,1:m+1,2)
end

function taylor_to_hermite1d3!(uhermite,utaylor,Hmat::Array{Float64,2},m::Int64)
  for j = 1:2*m+2
    @inbounds uhermite[j] = 0.0
    for i = 1:m+1
      @inbounds uhermite[j] += Hmat[j,i]*utaylor[i,1]+Hmat[j,m+1+i]*utaylor[i,2]
    end
  end
end

function taylor_to_hermite1d4!(uhermite,utaylor,Hmat::Array{Float64,2},m::Int64)
  for j = 1:2*m+2
    @inbounds uhermite[j] = 0.0
  end
  for i = 1:m+1
    for j = 1:2*m+2
      @inbounds uhermite[j] += Hmat[j,i]*utaylor[i,1]+Hmat[j,m+1+i]*utaylor[i,2]
    end
  end
end
