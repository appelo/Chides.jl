isdefined(@__MODULE__,:Chapter02) || include("common.jl")
using .Chapter02
function cardinal_basis(out;plots=true)
    r=cardinal_results()
    savecsv(joinpath(out,"cardinal_basis.csv"),["x";["B_left_$k" for k in 0:3];["B_right_$k" for k in 0:3]],[(r.x[i],r.curves[i,:]...) for i in eachindex(r.x)])
    if plots
        p=[]
        for side in 1:2
            push!(p,Chapter02.Plots.plot(r.x,r.curves[:,4side-3:4side],label=["k = $k" for _ in 1:1,k in 0:3],
                xlabel="z = (x - a) / h",ylabel="Basis value",title=side==1 ? "Left endpoint data" : "Right endpoint data"))
        end
        Chapter02.Plots.savefig(Chapter02.Plots.plot(p...;layout=(1,2)),joinpath(out,"hb_basis.pdf"))
    end
    r
end
if abspath(PROGRAM_FILE)==@__FILE__
    out=output_dir();loadplots();Base.invokelatest(cardinal_basis,out);Base.invokelatest(metadata,out)
end
