isdefined(@__MODULE__,:Chapter02) || include("common.jl")
using .Chapter02
function scaling_and_translation(out;plots=false)
    # Cubic midpoint formulas for scaled endpoint coefficients F0,F1,G0,G1.
    data=[1.2 -.7; .3 -.2]
    actual=reconstruct(data,1)
    F0,F1=data[:,1];G0,G1=data[:,2]
    formulas=[(F0+G0)/2+(F1-G1)/8,1.5*(G0-F0)-(F1+G1)/4,
              (G1-F1)/2,2*(F0-G0)+F1+G1]
    savecsv(joinpath(out,"cubic_midpoint.csv"),["coefficient","chides","formula","residual"],
        [(k,actual[k+1],formulas[k+1],abs(actual[k+1]-formulas[k+1])) for k in 0:3])
    println("Cubic midpoint formulas: maximum residual = ",maximum(abs.(actual.-formulas)))
    rows=scaling_results()
    savecsv(joinpath(out,"scaling_and_translation.csv"),["m","a","b","degree","center_fraction","maximum_error"],rows)
    println("Polynomial reproduction: maximum residual = ",maximum(row[end] for row in rows))
    translated=translation_results()
    savecsv(joinpath(out,"translation_and_rescaling.csv"),["old_center","old_width","new_center","new_width","coefficient_error","value_error"],translated)
    println("Translation and rescaling: maximum coefficient/value residual = ",maximum(max(r[5],r[6]) for r in translated))
    rows
end
if abspath(PROGRAM_FILE)==@__FILE__
    out=output_dir();scaling_and_translation(out);metadata(out)
end
