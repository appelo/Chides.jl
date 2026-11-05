module Chapter02
using Chides, TaylorSeries, LinearAlgebra, Printf, TOML
export hmap, reconstruct, evaluate, derivative, endpoint_data, exact_polynomial,
       polyderivative, polymultiply, polyintegral, inner, subtract, reference_step,
       reference_abs, barycentric, savecsv, output_dir, metadata, loadplots,
       cardinal_results, scaling_results, translation_results, step_results, absolute_results,
       smooth_results, projection_results, tensor_results

function hmap(m; a=-1.0, b=1.0, center=(a+b)/2, mode=0)
    H = zeros(2m+2,2m+2)
    hermitemap!(H,m,Float64(a),Float64(b),Float64(center),mode)
    H
end
function reconstruct(data, m; kwargs...)
    c = zeros(2m+2)
    taylor_to_hermite1d!(c,data,hmap(m;kwargs...),m)
    c
end
evaluate(c,x) = evalpoly(x,c)
function derivative(c,x,h=1.0)
    length(c)==1 && return zero(x)
    evalpoly(x,[(k-1)*c[k] for k in 2:length(c)])/h
end
function endpoint_data(f,a,b,m)
    h=b-a
    u=zeros(m+1,2)
    for (j,x) in enumerate((a,b))
        coefficients=f(x+Taylor1(Float64,m+1))
        for k in 0:m
            u[k+1,j]=h^k*coefficients[k]
        end
    end
    u
end
function polyderivative(c,q=1)
    q==0 && return copy(c)
    length(c)<=q && return [zero(eltype(c))]
    [c[k+1]*prod(big(k-q+1):big(k)) for k in q:length(c)-1]
end
function polymultiply(a,b)
    c=zeros(promote_type(eltype(a),eltype(b)),length(a)+length(b)-1)
    for i in eachindex(a), j in eachindex(b)
        c[i+j-1]+=a[i]*b[j]
    end
    c
end
polyintegral(c) = sum(2c[k+1]/(k+1) for k in 0:2:length(c)-1)
inner(a,b)=polyintegral(polymultiply(a,b))
function subtract(a,b)
    c=zeros(promote_type(eltype(a),eltype(b)),max(length(a),length(b)))
    c[1:length(a)].+=a; c[1:length(b)].-=b
    c
end
function exact_polynomial(f,m)
    R=Rational{BigInt}; A=zeros(R,2m+2,2m+2); rhs=zeros(R,2m+2)
    for (j,x) in enumerate((-1,1)), k in 0:m
        row=(j-1)*(m+1)+k+1
        for l in k:2m+1
            A[row,l+1]=binomial(big(l),k)*x^(l-k)
        end
        rhs[row]=evaluate(polyderivative(f,k),R(x))/factorial(big(k))
    end
    A\rhs
end
function reference_step(x,m)
    setprecision(256) do
        z=BigFloat(x)
        numerator=sum((-1)^k*BigFloat(binomial(big(m),k))*z^(2k+1)/(2k+1) for k in 0:m)
        denominator=sum((-1)^k*BigFloat(binomial(big(m),k))/(2k+1) for k in 0:m)
        -numerator/denominator
    end
end
function reference_abs(x,m)
    setprecision(256) do
        coefficient=BigFloat(1); value=coefficient; z=BigFloat(x)^2-1
        power=BigFloat(1)
        for k in 1:m
            coefficient*=(BigFloat(1)/2-k+1)/k
            power*=z; value+=coefficient*power
        end
        value
    end
end
function barycentric(x,nodes,values,weights)
    j=findfirst(==(x),nodes)
    j!==nothing && return values[j]
    terms=weights./(x.-nodes)
    sum(terms.*values)/sum(terms)
end
function savecsv(path,header,rows)
    open(path,"w") do io
        println(io,join(header,","))
        for row in rows
            println(io,join(map(x->x isa AbstractFloat ? @sprintf("%.17g",x) : string(x),row),","))
        end
    end
end
function output_dir(args=ARGS)
    index=findfirst(==("--output"),args)
    directory=index===nothing ? joinpath(@__DIR__,"output") : abspath(args[index+1])
    mkpath(directory); directory
end
function metadata(out)
    using_versions=Dict("Chides"=>string(pkgversion(Chides)),"TaylorSeries"=>string(pkgversion(TaylorSeries)))
    if isdefined(@__MODULE__,:Plots)
        using_versions["Plots"]=string(pkgversion(Plots))
    end
    repo=normpath(joinpath(@__DIR__,"../.."))
    commit=readchomp(`git -C $repo rev-parse HEAD`)
    dirty=!isempty(readchomp(`git -C $repo status --porcelain`))
    data=Dict("julia_version"=>string(VERSION),"chides_commit"=>commit,"working_tree_modified"=>dirty,
        "arithmetic"=>"Float64 reconstruction; 256-bit reference arithmetic; exact rational polynomial integration",
        "packages"=>using_versions,"parameters"=>Dict("basis_m"=>3,"step_m"=>[1,5,15],
        "absolute_m"=>[1,2,4,8,16],"absolute_reference_max_m"=>512,
        "smooth_m"=>[0,1,2,3],"smooth_cells"=>[2,4,8,16,32],"smooth_sample_points_per_cell"=>65,
        "projection_m"=>[0,1,2,3,4],"tensor_m"=>2,"tensor_cells"=>[8,6],"tensor_domain"=>[0,1,-1,1],"tensor_subintervals_per_cell"=>12,
        "basis_samples"=>501,"step_samples"=>4001,"absolute_samples"=>1001,
        "scaling_cells"=>[[-1.0,1.0],[0.3,0.6],[2.0,5.0]],"center_fractions"=>[0.0,0.27,0.5,1.0],
        "translation_centers"=>[-0.1,0.9,1.4],"translation_widths"=>[0.3,1.25,3.0]))
    open(joinpath(out,"metadata.toml"),"w") do io
        TOML.print(io,data;sorted=true)
    end
end
function loadplots()
    ENV["GKSwstype"]="100"
    @eval using Plots
    @eval begin
        Plots.gr()
        Plots.default(fontfamily="Computer Modern", linewidth=2, size=(700,350),
            guidefontsize=14,tickfontsize=12,legendfontsize=12,titlefontsize=16,margin=5Plots.mm)
    end
end

function cardinal_results()
    m=3; x=collect(range(0,1,length=501)); curves=zeros(length(x),2m+2)
    residual=0.0
    for j in 1:2m+2
        data=zeros(m+1,2); data[j]=1
        c=reconstruct(data,m)
        curves[:,j]=evaluate.(Ref(c),x.-0.5)
        for k in 0:m, (side,z) in enumerate((-.5,.5))
            got=evaluate(polyderivative(c,k),z)/factorial(k)
            residual=max(residual,abs(got-data[k+1,side]))
        end
    end
    (;x,curves,residual)
end
function scaling_results()
    rows=[]
    for m in 0:5, (a,b) in ((-1.0,1.0),(.3,.6),(2.0,5.0)), degree in 0:2m+1
        f(x)=((x-a)/(b-a))^degree
        data=endpoint_data(f,a,b,m)
        for fraction in (0.0,.27,.5,1.0)
            mode=fraction==0 ? -1 : fraction==1 ? 1 : 0
            center=a+fraction*(b-a)
            c=reconstruct(data,m;a,b,center,mode)
            error=maximum(abs(evaluate(c,(x-center)/(b-a))-f(x)) for x in range(a,b,length=41))
            push!(rows,(m,a,b,degree,fraction,error))
        end
    end
    rows
end
function translation_results()
    # Re-express the same degree-five polynomial using a new center and width.
    a,b=-0.4,1.6; center=(a+b)/2; h=b-a; m=2
    f(x)=1+2x-3x^2+x^3-.5x^4+.25x^5
    c=reconstruct(endpoint_data(f,a,b,m),m)
    rows=[]
    for newcenter in (-.1,.9,1.4), newwidth in (.3,1.25,3.0)
        sigma=(newcenter-center)/h
        translated=[(newwidth/h)^j*sum(binomial(k,j)*c[k+1]*sigma^(k-j) for k in j:2m+1) for j in 0:2m+1]
        reference=f(newcenter+Taylor1(Float64,2m+1))
        coefficient_error=maximum(abs(translated[j+1]-newwidth^j*reference[j]) for j in 0:2m+1)
        value_error=maximum(abs(evaluate(translated,(x-newcenter)/newwidth)-f(x)) for x in range(a,b,length=41))
        push!(rows,(center,h,newcenter,newwidth,coefficient_error,value_error))
    end
    rows
end
function step_results()
    x=collect(range(-1,1,length=4001)); curves=[]; rows=[]
    for m in (1,5,15)
        data=zeros(m+1,2);data[1,:]=[1,-1]
        c=reconstruct(data,m); values=evaluate.(Ref(c),x./2)
        n=2m+1; nodes=[-1+2BigFloat(k)/n for k in 0:n]; data_l=-sign.(nodes)
        weights=[BigFloat((-1)^k*binomial(big(n),k)) for k in 0:n]
        lagrange=[Float64(setprecision(256) do
            barycentric(BigFloat(z),BigFloat.(nodes),BigFloat.(data_l),BigFloat.(weights))
        end) for z in x]
        discrepancy=maximum(abs(values[i]-Float64(reference_step(x[i],m))) for i in 1:20:length(x))
        endpoint_error=max(abs(values[1]-1),abs(values[end]+1))
        variation=sum(abs.(diff(values)))
        overshoot=max(0.0,maximum(abs.(values))-1)
        push!(rows,(m,endpoint_error,variation,overshoot,discrepancy,maximum(abs.(lagrange))))
        push!(curves,(m,values,lagrange))
    end
    (;x,curves,rows)
end
function absolute_results()
    x=collect(range(-1,1,length=1001)); curves=[];rows=[]
    for m in (1,2,4,8,16)
        data=zeros(m+1,2);data[1,:].=1;data[2,:]=[-2,2]
        c=reconstruct(data,m); values=evaluate.(Ref(c),x./2)
        exact=BigFloat(binomial(big(2m),m))/BigFloat(4)^m
        discrepancy=maximum(abs(values[i]-Float64(reference_abs(x[i],m))) for i in 1:10:length(x))
        push!(rows,(m,c[1],Float64(exact),discrepancy));push!(curves,(m,values))
    end
    reference=[(m,Float64(BigFloat(binomial(big(2m),m))/BigFloat(4)^m)) for m in 1:512]
    (;x,curves,rows,reference)
end
function smooth_results()
    rows=[]
    for m in 0:3, N in (2,4,8,16,32)
        h=1/N; errors=zeros(2)
        for cell in 0:N-1
            a=cell*h;b=a+h;center=(a+b)/2
            c=reconstruct(endpoint_data(x->sin(2pi*x),a,b,m),m)
            for z in range(-.5,.5,length=65)
                x=center+h*z
                errors[1]=max(errors[1],abs(evaluate(c,z)-sin(2pi*x)))
                errors[2]=max(errors[2],abs(derivative(c,z,h)-2pi*cos(2pi*x)))
            end
        end
        push!(rows,(m,N,h,errors[1],errors[2],errors[1]<1e-12,errors[2]<1e-11))
    end
    rows
end
function projection_results()
    rows=[]
    for m in 0:4
        R=Rational{BigInt}; f=zeros(R,2m+4);f[2]=1;f[2m+3]=1;f[2m+4]=1
        exact=exact_polynomial(f,m)
        data=zeros(m+1,2)
        for (side,x) in enumerate((-1,1)), k in 0:m
            data[k+1,side]=Float64(2^k*evaluate(polyderivative(f,k),R(x))/factorial(big(k)))
        end
        c=reconstruct(data,m)
        p=[BigFloat(c[k+1])/BigFloat(2)^k for k in 0:2m+1]
        df=BigFloat.(polyderivative(f,m+1)); dp=polyderivative(p,m+1); de=subtract(df,dp)
        ef=inner(df,df); ep=inner(dp,dp); ee=inner(de,de)
        orth=abs(inner(de,dp))/sqrt(ef*max(ep,eps(BigFloat)))
        pyth=abs(ef-ep-ee)/ef
        w=R[1]; for k in 1:m+1;w=polymultiply(w,R[-1,0,1]);end
        dw=BigFloat.(polyderivative(w,m+1)); pert=subtract(dp,-dw)
        minimal=abs(inner(pert,pert)-ep-inner(dw,dw))/max(ef,inner(dw,dw))
        # The independent exact calculation must satisfy the identities identically.
        dref=polyderivative(exact,m+1); dexact=polyderivative(f,m+1)
        @assert inner(subtract(dexact,dref),dref)==0
        @assert inner(dexact,dexact)==inner(dref,dref)+inner(subtract(dexact,dref),subtract(dexact,dref))
        @assert inner(dref,polyderivative(w,m+1))==0
        coefficient_error=maximum(abs.(subtract(p,BigFloat.(exact))))/max(1,maximum(abs.(exact)))
        push!(rows,(m,Float64(orth),Float64(pyth),Float64(minimal),Float64(coefficient_error)))
    end
    rows
end
function tensor_results(;polynomial=false, omit_mixed=false)
    m=2;nx=9;ny=7;nref=12;hx=1/(nx-1);hy=2/(ny-1)
    xp=range(0,1,length=nx);yp=range(-1,1,length=ny)
    up=zeros(m+1,m+1,nx,ny)
    for j in 1:ny,i in 1:nx,k in 0:m,l in 0:m
        if polynomial
            up[k+1,l+1,i,j]=hx^k*hy^l*binomial(5,k)*binomial(4,l)*xp[i]^(5-k)*yp[j]^(4-l)
        else
            up[k+1,l+1,i,j]=hx^k*hy^l*exp(xp[i])*pi^l*cos(pi*yp[j]+l*pi/2)/(factorial(k)*factorial(l))
        end
        omit_mixed && k>0 && l>0 && (up[k+1,l+1,i,j]=0)
    end
    values=zeros((nx-1)*nref+1,(ny-1)*nref+1)
    primal2plot2d!(values,up,hmap(m),nref)
    x=collect(range(0,1,length=size(values,1)));y=collect(range(-1,1,length=size(values,2)))
    truth=polynomial ? [a^5*b^4 for a in x,b in y] : [exp(a)*cos(pi*b) for a in x,b in y]
    (;x,y,values,error=abs.(values.-truth),maximum_error=maximum(abs.(values.-truth)))
end
end
