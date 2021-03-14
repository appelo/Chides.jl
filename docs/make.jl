push!(LOAD_PATH,"../src/")
using Documenter, Example, Chides

makedocs(sitename="Chides Documentation")

makedocs(
    modules = [Chides],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    #format = Documenter.LaTeX(platform = "docker"), # errors with no such file or directory???
    #clean = false,
    sitename="Chides.jl",
    authors = "Daniel Appelo, Thomas Hagstrom, and contributors.",
    pages = [
        "Home" => "index.md",
        "Workflow" => "workflow.md",
        "Examples" => "examples.md",
        "Types" => "types.md",
        #"Methods" => "methods.md",
        "Index" => "function-index.md",
    ],
    doctest = false, # set to true once initial implementation is done
    source = "src",
)

#deploydocs(
#    repo = "github.com/appelo/Chides.jl",
#    target = "build",
#    devurl="docs",
#)
