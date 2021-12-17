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
using Makie, CairoMakie, CSV, DataFrames, SnpArrays, Statistics, ColorSchemes

const GM = GeneticsMakie
isdir("data") || mkdir("data")
isdir("figs") || mkdir("figs")

# Download the latest GENCODE annotation
gencode = let
    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/gencode.v39lift37.annotation.gtf.gz"
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
range1 = start - 1e6
range2 = stop + 1e6

# Visualize 1 Mb window around KMT2E
set_theme!(font = "Arial")
let
    f = Figure(resolution = (306, 250))
    ax = Axis(f[1, 1])
    rs = GM.plotgenes!(ax, chr, range1, range2, gencode; height = 0.1)
    vspan!(ax, start, stop, color = (:gray, 0.2))
    Label(f[1, 1, Bottom()], "~$(round(range1 / 1e6; digits = 1)) Mb", textsize = 6, halign = :left)
    Label(f[1, 1, Bottom()], "Chr $(chr)", textsize = 6, halign = :center)
    Label(f[1, 1, Bottom()], "~$(round(range2 / 1e6; digits = 1)) Mb", textsize = 6, halign = :right)
    rowsize!(f.layout, 1, rs)
    save("figs/KMT2E-gene.png", f, px_per_unit = 4)
end
```
<p align="center"><img width="50%" style="border-radius: 5px;" src="figs/KMT2E-gene.png"></p>

```julia
# Visualize KMT2E isoforms
let
    f = Figure(resolution = (306, 306))
    ax = Axis(f[1, 1])
    rs, chr, range1, range2 = GM.plotisoforms!(ax, "KMT2E", gencode; height = 0.1)
    Label(f[1, 1, Bottom()], "~$(round(range1 / 1e6; digits = 1)) Mb", textsize = 6, halign = :left)
    Label(f[1, 1, Bottom()], "Chr $(chr)", textsize = 6, halign = :center)
    Label(f[1, 1, Bottom()], "~$(round(range2 / 1e6; digits = 1)) Mb", textsize = 6, halign = :right)
    rowsize!(f.layout, 1, rs)
    save("figs/KMT2E-isoform.png", f, px_per_unit = 4)
end
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
    isfile("data/$(meta)") || download(url, "data/$(meta)")
    # Subset data to the genomic region of interest and European samples
    kgp = SnpData("data/$(replace(vcf, ".vcf.gz" => ""))")
    meta = CSV.read("data/$(meta)", DataFrame)
    eur = meta.sample[meta.super_pop .== "EUR"]
    colinds = findall((kgp.snp_info.position .>= range1) .& (kgp.snp_info.position .<= range2))
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
size(kgp) # 503 European individuals + 5,574 SNPs

geno = convert(Matrix{Float64}, kgp.snparray)
LD = cor(geno, dims = 1)
LD = LD.^2 

# Visualize LD for KMT2E locus
@time let
    f = Figure(resolution = (306, 200))
    ax = Axis(f[1, 1])
    GM.plotld!(ax, LD)
    Label(f[1, 1, Top()], "~$(round(range1 / 1e6; digits = 1)) Mb", textsize = 6, halign = :left)
    Label(f[1, 1, Top()], "Chr $(chr)", textsize = 6, halign = :center)
    Label(f[1, 1, Top()], "~$(round(range2 / 1e6; digits = 1)) Mb", textsize = 6, halign = :right)
    rowsize!(f.layout, 1, 137)
    Colorbar(f[2, 1], limits = (0, 1), tellwidth = false, ticklabelsize = 6,
        colormap = cgrad(ColorSchemes.Blues_9, 9, categorical = true),
        label = "LD", labelsize = 6, vertical = false, flipaxis = false,
        width = 40, spinewidth = 0.5, tickwidth = 0.5, height = 5, ticksize = 3)
    rowgap!(f.layout, 5)
    save("figs/KMT2E-LD.png", f, px_per_unit = 4)
end
```
<p align="center"><img width="60%" style="border-radius: 5px;" src="figs/KMT2E-LD.png"></p>

```julia
# Download GWAS summary statistics for psychiatric disorders
gwas = let
    gwas = Dict("scz" => ("https://figshare.com/ndownloader/files/28169757", "PGC3_SCZ_wave3_public.v2.tsv.gz"),
        "bd" => ("https://figshare.com/ndownloader/files/26603681", "pgc-bip2021-all.vcf.tsv.gz"),
        "asd" => ("https://figshare.com/ndownloader/files/28169292", "iPSYCH-PGC_ASD_Nov2017.gz"))
    for key in keys(gwas)
        url = gwas[key][1]
        file = gwas[key][2]
        isfile("data/$(file)") || download(url, "data/$(file)")
    end
    scz = CSV.File("data/$(gwas["scz"][2])"; delim = "\t") |> DataFrame
    dropmissing!(scz) # row 5,328,757 is funky and contains missing values
    bd = CSV.File("data/$(gwas["bd"][2])", delim = "\t", skipto = 74, 
        header = ["CHROM", "POS", "ID", "A1", "A2", "BETA", "SE", "P", "NGT", 
            "FCAS", "FCON", "IMPINFO", "NEFFDIV2", "NCAS", "NCON", "DIRE"]) |> DataFrame
    asd = CSV.File("data/$(gwas["asd"][2])"; delim = "\t") |> DataFrame
    [scz, bd, asd]
end
titles = ["Schizophrenia (PGC3)", "Bipolar (Mullins et al. 2021)", "Autism (Grove et al. 2019)"]

# Visualize GWAS results for KMT2E locus
let
    n = length(gwas)
    f = Figure(resolution = (306, 350))
    axs = [Axis(f[i, 1]) for i in 1:(n + 1)]
    for i in 1:n
        GM.plotlocus!(axs[i], chr, range1, range2, gwas[i]; colorld = true, ref = kgp, ymax = 18)
        rowsize!(f.layout, i, 30)
        Label(f[i, 1, Top()], "$(titles[i])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
    end
    rs = GM.plotgenes!(axs[n + 1], chr, range1, range2, gencode; height = 0.1)
    rowsize!(f.layout, n + 1, rs)
    Label(f[n + 1, 1, Bottom()], "~$(round(range1 / 1e6; digits = 1)) Mb", textsize = 6, halign = :left)
    Label(f[n + 1, 1, Bottom()], "Chr $(chr)", textsize = 6, halign = :center)
    Label(f[n + 1, 1, Bottom()], "~$(round(range2 / 1e6; digits = 1)) Mb", textsize = 6, halign = :right)
    Colorbar(f[1:n, 2], limits = (0, 1), ticks = 0:1:1, height = 20,
        colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
        tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
        labelsize = 6, width = 5, spinewidth = 0.5)
    Label(f[1:n, 0], text = "-log[p]", textsize = 6, rotation = pi / 2)
    rowgap!(f.layout, 5)
    colgap!(f.layout, 5)
    for i in 1:(n + 1)
        vlines!(axs[i], start, color = (:gold, 0.9), linewidth = 0.5)
    end
    save("figs/KMT2E-locuszoom.png", f, px_per_unit = 4)
end
```
<p align="center"><img width="70%" style="border-radius: 5px;" src="figs/KMT2E-locuszoom.png"></p>

```julia
# Visualize QQ plot of P values
let
    f = Figure(resolution = (612, 255))
    axs = [Axis(f[2, i]) for i in 1:3]
    for i in eachindex(titles)
        GM.plotqq!(axs[i], gwas[i]; xlabel = "", ylabel = "", ystep = 5)
        ylims!(axs[i], 0, 40)
        i > 1 ? hideydecorations!(axs[i]) : nothing
    end
    for (i, title) in enumerate(titles)
        Box(f[1, i], color = :gray90)
        Label(f[1, i], title, tellwidth = false, textsize = 8, padding = (3, 0, 3, 3))
    end
    Label(f[3, 1:length(titles)], text = "Expected -log[p]", textsize = 8)
    Label(f[2, 0], text = "Observed -log[p]", textsize = 8, rotation = pi / 2, tellheight = false)
    rowsize!(f.layout, 2, Aspect(2, 1))
    colgap!(f.layout, 5)
    rowgap!(f.layout, 1, 0)
    rowgap!(f.layout, 2, 5)
    save("figs/QQ.png", f, px_per_unit = 4)
end
```
<p align="center"><img width="65%" style="border-radius: 5px;" src="figs/QQ.png"></p>