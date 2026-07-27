using Documenter, Chides

DocMeta.setdocmeta!(Chides, :DocTestSetup, :(using Chides); recursive = true)

makedocs(
    modules = [Chides],
    sitename = "Chides.jl",
    authors = "Daniel Appelo, Thomas Hagstrom, and contributors.",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://appelo.github.io/Chides.jl/stable/",
        edit_link = "master",
    ),
    pages = [
        "Home" => "index.md",
        "Workflow" => "workflow.md",
        "Examples" => "examples.md",
        "Index" => "function-index.md",
    ],
)

deploydocs(
    repo = "github.com/appelo/Chides.jl",
    devbranch = "master",
)
