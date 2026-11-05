isdefined(@__MODULE__,:Chapter02) || include("common.jl")
using .Chapter02
function tensor_interpolation(out;plots=true)
    r=tensor_results(); omitted=tensor_results(omit_mixed=true)
    savecsv(joinpath(out,"tensor_interpolation.csv"),["x","y","hermite","absolute_error"],[(r.x[i],r.y[j],r.values[i,j],r.error[i,j]) for j in eachindex(r.y) for i in eachindex(r.x)])
    savecsv(joinpath(out,"tensor_diagnostics.csv"),["complete_data_error","omitted_mixed_error"],[(r.maximum_error,omitted.maximum_error)])
    if plots
        p=Chapter02.Plots.heatmap(r.x,r.y,r.values';xlabel="x",ylabel="y",title="Hermite reconstruction",color=:viridis)
        q=Chapter02.Plots.heatmap(r.x,r.y,(1e5.*r.error)';xlabel="x",ylabel="y",title="Absolute error / 1e-5",color=:inferno)
        Chapter02.Plots.savefig(Chapter02.Plots.plot(p,q;layout=(1,2),size=(750,350)),joinpath(out,"hb_tensor.pdf"))
    end
    r
end
if abspath(PROGRAM_FILE)==@__FILE__
    out=output_dir();loadplots();Base.invokelatest(tensor_interpolation,out);Base.invokelatest(metadata,out)
end
