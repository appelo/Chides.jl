isdefined(@__MODULE__,:Chapter02) || include("common.jl")
using .Chapter02
function step_function(out;plots=true)
    r=step_results()
    savecsv(joinpath(out,"step_diagnostics.csv"),["m","endpoint_error","sampled_total_variation","overshoot","reference_discrepancy","max_abs_lagrange"],r.rows)
    savecsv(joinpath(out,"step_curves.csv"),["m","x","hermite","lagrange"],[(m,r.x[i],v[i],l[i]) for (m,v,l) in r.curves for i in eachindex(r.x)])
    if plots
        panels=[]
        for (m,v,l) in r.curves
            p=Chapter02.Plots.plot(r.x,l;label="Lagrange",xlabel="x",ylabel="Value",title="m = $m: full range",color=2,xticks=[-1,0,1],legendfontsize=10,titlefontsize=14)
            Chapter02.Plots.plot!(p,r.x,v;label="Hermite",color=1)
            if m==15
                Chapter02.Plots.plot!(p;yticks=([-1e6,0,1e6],["-1e6","0","1e6"]))
            end
            push!(panels,p)
        end
        for (m,v,l) in r.curves
            p=Chapter02.Plots.plot(r.x,l;label="Lagrange",xlabel="x",ylabel="Value",title="m = $m: jump detail",color=2,xlims=(-.5,.5),ylims=(-1.5,1.5),xticks=[-.5,0,.5],yticks=[-1,0,1],legendfontsize=10,titlefontsize=14)
            Chapter02.Plots.plot!(p,r.x,v;label="Hermite",color=1)
            Chapter02.Plots.plot!(p,r.x,-sign.(r.x);label="Step",color=:black,linestyle=:dash)
            push!(panels,p)
        end
        Chapter02.Plots.savefig(Chapter02.Plots.plot(panels[[1,4,2,5,3,6]]...;layout=(3,2),size=(700,800)),joinpath(out,"hb_step.pdf"))
    end
    r
end
if abspath(PROGRAM_FILE)==@__FILE__
    out=output_dir();loadplots();Base.invokelatest(step_function,out);Base.invokelatest(metadata,out)
end
