include("common.jl")
using .Chapter02
for name in ("cardinal_basis","scaling_and_translation","step_function","absolute_value","smooth_convergence","projection_identity","tensor_interpolation")
    include(name*".jl")
end
out=output_dir();loadplots()
cardinal_basis(out);scaling_and_translation(out);step_function(out);absolute_value(out)
smooth_convergence(out);projection_identity(out);tensor_interpolation(out);metadata(out)
index=findfirst(==("--book-root"),ARGS)
if index!==nothing
    index<length(ARGS) || error("--book-root requires a directory")
    book=abspath(ARGS[index+1]);isfile(joinpath(book,"HermiteBook.tex")) || error("Invalid book root")
    destination=joinpath(book,"Chapter_HB_definition","figs")
    for name in ("hb_basis.pdf","hb_step.pdf","hb_abs.pdf","hb_convergence.pdf","hb_tensor.pdf","hb_projection_table.tex")
        cp(joinpath(out,name),joinpath(destination,name);force=true)
    end
    println("Exported six Chapter 2 assets to ",destination)
end
println("Chapter 2 outputs: ",out)
