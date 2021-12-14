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

```julia
kgp = let
    # Download 1000 Genomes data for a single chromosome
    beagle = "http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a"
    url = joinpath(beagle, "b37.vcf/chr$(chr).1kg.phase3.v5a.vcf.gz")
    vcf = last(split(url, "/"))
    isfile("data/$(vcf)") || download(url, "data/$(vcf)")
    # Convert vcf file to plink bed file (this step takes a while)
    @time isfile("data/$(replace(vcf, ".vcf.gz" => ".bed"))") || vcf2plink("data/$(vcf)", "data/$(replace(vcf, ".vcf.gz" => ""))")
    # Download sample metadata
    url = joinpath(beagle, "/sample_info/integrated_call_samples_v3.20130502.ALL.panel")
    meta = last(split(url, "/")) 
    isfile("data/$(meta)") || download(url2, "data/$(meta)")
    # Subset data to the genomic region of interest and European samples
    kgp = SnpData("data/$(replace(vcf, ".vcf.gz" => ""))")
    meta = CSV.read("data/$(meta)", DataFrame)
    eur = meta.sample[meta.super_pop .== "EUR"]
    colinds = findall((kgp.snp_info.position .>= start - 1e6) .& (kgp.snp_info.position .<= stop + 1e6))
    rowinds = findall(in(eur), kgp.person_info.iid)
    file = replace(vcf, ".vcf.gz" => ".eur")
    SnpArrays.filter(kgp, rowinds, colinds; des = "data/$(file)")
    # Apply minor allele frequency > 0.05 filter
    kgp = SnpData("data/$(file)")
    colinds = SnpArrays.filter(kgp.snparray; min_maf = 0.05)[2]
    file = file * ".maf0.05"
    SnpArrays.filter(kgp, trues(size(kgp)[1]), colinds; des = "data/$(file)")
    SnpData("data/$(file)")
end 
size(kgp) # 503 European individuals + 5,346 SNPs

geno = convert(Matrix{Float64}, kgp.snparray)
LD = cor(geno, dims = 1)
LD = LD.^2 

# Visualize LD for KMT2E locus
@time GM.plotld(LD, filename = "figs/KMT2E-LD"; 
    xlabel = (kgp.snp_info.chromosome[1], kgp.snp_info.position[1], kgp.snp_info.position[end]))
```
<p align="center"><img width="60%" style="border-radius: 5px;" src="figs/KMT2E-LD.png"></p>
