# GeneticsMakie

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mmkim1210.github.io/GeneticsMakie.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmkim1210.github.io/GeneticsMakie.jl/dev)
[![Build Status](https://github.com/mmkim1210/GeneticsMakie.jl/workflows/CI/badge.svg)](https://github.com/mmkim1210/GeneticsMakie.jl/actions)
[![Coverage](https://codecov.io/gh/mmkim1210/GeneticsMakie.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmkim1210/GeneticsMakie.jl)

## Installation
```julia
julia>]
pkg> add https://github.com/mmkim1210/GeneticsMakie.jl.git
```

## Examples
```julia
using GeneticsMakie
using CSV, DataFrames, DataFramesMeta, Chain, SnpArrays, Statistics

const GM = GeneticsMakie
isdir("data") || mkdir("data")
isdir("figs") || mkdir("figs")

# Download the latest GENCODE annotation
gencode = let 
    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gtf.gz"
    file = last(split(url, "/"))
    isfile("data/$(file)") || download(url, "data/$(file)")
    header = ["seqnames", "source", "feature", "start", "end", "score", "strand", "phase", "info"]
    CSV.File("data/$(file)"; delim = "\t", skipto = 6, header = header) |> DataFrame
end

# Parse GENCODE 
GM.parsegtf!(gencode)

# Focus on KMT2E gene as an example
chr, start, stop = let 
    ind = findfirst(isequal("KMT2E"), gencode.gene_name)
    gencode.seqnames[ind], gencode.start[ind], gencode[ind, :end]
end

# Visualize 1 Mb window around KMT2E
GM.plotgenes(chr, start - 1e6, stop + 1e6, gencode; filename = "figs/KMT2E-gene", height = 0.15)
```
<p align="center"><img width="50%" style="border-radius: 5px;" src="figs/KMT2E-gene.png"></p>

```julia
# Visualize KMT2E isoforms
GM.plotisoforms("KMT2E", gencode; filename = "figs/KMT2E-isoform", height = 0.1)
```
<p align="center"><img width="50%" style="border-radius: 5px;" src="figs/KMT2E-isoform.png"></p>