

function test_hmap(m)
    Hmat = zeros(Float64,2*(m+1),2*(m+1))
    xl = -1.0
    xr = 1.0
    xc = 0.0
    icase = 0
    hermitemap!(Hmat,m,xl,xr,xc,icase)
    return isapprox(Hmat,[[0.5,-1.0] [0.5,1.0]])
    # return true
end
