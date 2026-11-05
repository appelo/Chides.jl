isdefined(@__MODULE__,:Chapter02) || include("common.jl")
using .Chapter02
function smooth_convergence(out;plots=true)
    rows=smooth_results()
    savecsv(joinpath(out,"smooth_convergence.csv"),["m","cells","h","value_error","derivative_error","value_roundoff","derivative_roundoff"],rows)
    rates=[]
    for m in 0:3,kind in 1:2
        resolved=filter(r->r[1]==m && !r[5+kind],rows)
        a,b=resolved[end-1:end]
        push!(rates,(m,kind==1 ? "value" : "derivative",a[2],b[2],log2(a[3+kind]/b[3+kind]),2m+3-kind))
    end
    savecsv(joinpath(out,"smooth_rates.csv"),["m","quantity","coarse_cells","fine_cells","observed_order","expected_order"],rates)
    if plots
        panels=[]
        for kind in 1:2
            p=Chapter02.Plots.plot(;xscale=:log10,yscale=:log10,xlabel="Number of cells",xticks=([2,4,8,16,32],string.([2,4,8,16,32])),ylabel="Maximum sampled error",title=kind==1 ? "Function value" : "First derivative",legend=:outerbottom,legend_column=2,legendfontsize=10)
            for m in 0:3
                r=filter(r->r[1]==m,rows)
                Chapter02.Plots.plot!(p,getindex.(r,2),getindex.(r,3+kind);marker=:circle,label="m=$m, order $(2m+3-kind)")
            end
            push!(panels,p)
        end
        Chapter02.Plots.savefig(Chapter02.Plots.plot(panels...;layout=(1,2),size=(700,430)),joinpath(out,"hb_convergence.pdf"))
    end
    (;rows,rates)
end
if abspath(PROGRAM_FILE)==@__FILE__
    out=output_dir();loadplots();Base.invokelatest(smooth_convergence,out);Base.invokelatest(metadata,out)
end
