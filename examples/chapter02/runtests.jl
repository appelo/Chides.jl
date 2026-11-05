include("common.jl")
using .Chapter02, Chides, Test, LinearAlgebra
@testset "Chapter 2 numerical checks" begin
    @test maximum(r[end] for r in scaling_results())<1e-10
    @test maximum(max(r[5],r[6]) for r in translation_results())<1e-10
    @test cardinal_results().residual<1e-11
    # Explicit endpoint cardinal formulas, independent of the Fortran map.
    m=3
    basis(z,k)=z^k*(1-z)^(m+1)*sum(binomial(m+j,j)*z^j for j in 0:m-k)
    r=cardinal_results()
    for k in 0:m
        @test maximum(abs.(r.curves[:,k+1].-[basis(z,k) for z in r.x]))<1e-12
        @test maximum(abs.(r.curves[:,m+k+2].-[(-1)^k*basis(1-z,k) for z in r.x]))<1e-12
    end
    data=[1.2 -.7; .3 -.2]
    c=reconstruct(data,1)
    @test c[1]≈(data[1,1]+data[1,2])/2+(data[2,1]-data[2,2])/8
    @test c[2]≈1.5*(data[1,2]-data[1,1])-(data[2,1]+data[2,2])/4
    @test c[3]≈(data[2,2]-data[2,1])/2
    @test c[4]≈2*(data[1,1]-data[1,2])+data[2,1]+data[2,2]
    for m in 0:5
        @test hmap(m;mode=-1)[1:m+1,1:m+1]≈Matrix{Float64}(I,m+1,m+1)
        @test hmap(m;mode=1)[1:m+1,m+2:2m+2]≈Matrix{Float64}(I,m+1,m+1)
    end
    for r in step_results().rows
        @test r[2]<1e-9
        @test abs(r[3]-2)<1e-9
        @test r[4]<1e-9
        @test r[5]<1e-9
    end
    for r in absolute_results().rows
        @test abs(r[2]-r[3])<1e-9
        @test r[4]<1e-9
    end
    nodes=BigFloat[-1,0,1];values=nodes.^2;weights=BigFloat[1,-2,1]
    @test all(barycentric(nodes[i],nodes,values,weights)==values[i] for i in eachindex(nodes))
    @test abs(barycentric(BigFloat("0.27"),nodes,values,weights)-BigFloat("0.27")^2)<big"1e-70"
    rows=smooth_results()
    for m in 0:3,kind in 1:2
        resolved=filter(r->r[1]==m && !r[5+kind],rows)
        a,b=resolved[end-1:end]
        @test abs(log2(a[3+kind]/b[3+kind])-(2m+3-kind))<.25
    end
    @test maximum(maximum(r[2:end]) for r in projection_results())<1e-10
    @test tensor_results(polynomial=true).maximum_error<1e-12
    @test tensor_results(polynomial=true,omit_mixed=true).maximum_error>1e-5
    @test tensor_results().maximum_error<1e-4
    @test tensor_results(omit_mixed=true).maximum_error>10tensor_results().maximum_error
    m=3;nx=5;nref=8;h=1/(nx-1);x=range(0,1,length=nx)
    up=zeros(m+1,nx)
    for i in 1:nx,k in 0:m
        up[k+1,i]=h^k*binomial(7,k)*x[i]^(7-k)
    end
    values=zeros((nx-1)*nref+1)
    Chides.primal2plot1d!(values,up,hmap(m),nref)
    @test maximum(abs.(values.-collect(range(0,1,length=length(values))).^7))<1e-12
end
