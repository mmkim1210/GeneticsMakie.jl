# Munging summary statistics
GWAS summary statistics come in a variety of shapes and flavors, so harmonizing them
is crucial in making our lives easier when trying to visualize their results. We can 
take a peak at GWAS results for height and weight, the two classic anthropometric traits. 

```julia
using Pkg
Pkg.add(["GeneticsMakie", "CSV", "DataFrames"])
```

```julia
using GeneticsMakie, CSV, DataFrames
gwas = Dict(
    "height" => "https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz",
    "weight" => "https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz"
)
isdir("data/gwas") || mkdir("data/gwas")
df = DataFrame[]
for key in keys(gwas)
    url = gwas[key]
    isfile("data/gwas/$(basename(url))") || download(url, "data/gwas/$(basename(url))")
    push!(df, CSV.read("data/gwas/$(basename(url))", DataFrame; comment = "##", missingstring = ["NA"]))
end
```

To harmonize summary statistics, we run a single command `GeneticsMakie.mungesumstats!`. 
For all downstream plotting functions, we note that GWAS summary statistics should be pre-harmonized.
```julia
GeneticsMakie.mungesumstats!(df)
```
```@docs
mungesumstats!
```

To find GWAS loci for each phenotype that are separated by at least 1 Mb, 
we can execute the following command.
```julia
GeneticsMakie.findgwasloci(df[1])
GeneticsMakie.findgwasloci(df[2])
```

To find every GWAS loci across multiple phenotypes that are separated by at least 1 Mb,
we can instead run the following command.
```julia
GeneticsMakie.findgwasloci(df)
```
```@docs
findgwasloci
```

Such an exhaustive list of loci can then be iterated through and visualized by 
[Plotting LocusZooom](@ref).

We can also find the closest gene to each index SNP in GWAS loci.
```julia
loci = GeneticsMakie.findgwasloci(df[1])
GeneticsMakie.findclosestgene(loci, gencode)
GeneticsMakie.findclosestgene(loci, gencode; start = true) # closest gene from gene start site
GeneticsMakie.findclosestgene(loci, gencode; proteincoding = true) # closest "protein-coding" gene
```
```@docs
findclosestgene
```