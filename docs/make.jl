using GeneticsMakie
using Documenter
# using other stuff

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
        sidebar_sitename = false,
    ),
    pages=[
        "Home" => "index.md",
        "Examples" =>
            ["Parsing GENCODE" => "examples/gtf.md",
            "Munging summary statistics" => "examples/summary.md",
            "Plotting genes" => "examples/genes.md",
            "Plotting isoforms" => "examples/isoforms.md",
            "Plotting LocusZoom" => "examples/locus.md",
            "Plotting GWAS" => "examples/gwas.md",
            "Plotting TWAS" => "examples/twas.md",
            "Plotting ChIP-Seq & ATAC-Seq" => "examples/chip.md"],
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/mmkim1210/GeneticsMakie.jl",
)
