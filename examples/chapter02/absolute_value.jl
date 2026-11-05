isdefined(@__MODULE__,:Chapter02) || include("common.jl")
using .Chapter02
function absolute_value(out;plots=true)
    r=absolute_results()
    savecsv(joinpath(out,"absolute_diagnostics.csv"),["m","chides_cusp","exact_cusp","reference_discrepancy"],r.rows)
    savecsv(joinpath(out,"absolute_reference.csv"),["m","exact_cusp"],r.reference)
    savecsv(joinpath(out,"absolute_curves.csv"),["m","x","hermite"],[(m,r.x[i],v[i]) for (m,v) in r.curves for i in eachindex(r.x)])
    if plots
        p=Chapter02.Plots.plot(r.x,abs.(r.x);label="|x|",color=:black,linestyle=:dash,xlabel="x",ylabel="Value",legend=:outerbottom,legend_column=2,legendfontsize=10)
        for (m,v) in r.curves
            Chapter02.Plots.plot!(p,r.x,v;label="m = $m")
        end
        q=Chapter02.Plots.plot(first.(r.reference),last.(r.reference);label="Exact formula",xscale=:log10,yscale=:log10,xlabel="m",ylabel="Error at x = 0",legend=:outerbottom,legendfontsize=10)
        Chapter02.Plots.plot!(q,1:512,1 ./sqrt.(pi.*(1:512));label="1 / √(π m)",linestyle=:dash)
        Chapter02.Plots.scatter!(q,first.(r.rows),getindex.(r.rows,2);label="Chides (m <= 16)",markersize=4)
        Chapter02.Plots.savefig(Chapter02.Plots.plot(p,q;layout=(1,2),size=(700,440)),joinpath(out,"hb_abs.pdf"))
    end
    r
end
if abspath(PROGRAM_FILE)==@__FILE__
    out=output_dir();loadplots();Base.invokelatest(absolute_value,out);Base.invokelatest(metadata,out)
end
