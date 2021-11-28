using GeneticsMakie
using Documenter

DocMeta.setdocmeta!(GeneticsMakie, :DocTestSetup, :(using GeneticsMakie); recursive=true)

makedocs(;
    modules=[GeneticsMakie],
    authors="Minsoo <mmkim1210@gmail.com> and contributors",
    repo="https://github.com/mmkim1210/GeneticsMakie.jl/blob/{commit}{path}#{line}",
    sitename="GeneticsMakie.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mmkim1210.github.io/GeneticsMakie.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mmkim1210/GeneticsMakie.jl",
)
