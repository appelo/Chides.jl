isdefined(@__MODULE__,:Chapter02) || include("common.jl")
using .Chapter02, Printf
function projection_identity(out;plots=false)
    rows=projection_results()
    savecsv(joinpath(out,"projection_identity.csv"),["m","orthogonality","pythagorean","minimality","coefficient_error"],rows)
    open(joinpath(out,"hb_projection_table.tex"),"w") do io
        println(io,raw"\begin{tabular}{rcccc}")
        println(io,raw"\toprule")
        println(io,raw"$m$ & $R_{\rm orth}$ & $R_{\rm Pyth}$ & $R_{\rm min}$ & Coefficients ",repeat("\\",2))
        println(io,raw"\midrule")
        for row in rows
            @printf(io,"%d & %.1e & %.1e & %.1e & %.1e \\\\\n",row...)
        end
        println(io,raw"\bottomrule")
        println(io,raw"\end{tabular}")
    end
    rows
end
if abspath(PROGRAM_FILE)==@__FILE__
    out=output_dir();projection_identity(out);metadata(out)
end
