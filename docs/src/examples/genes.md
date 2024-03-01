# Plotting genes

After [Parsing GENCODE](@ref), we can start plotting gene bodies. 
[__GeneticsMakie.jl__](https://github.com/mmkim1210/GeneticsMakie.jl) is transparent 
in that it shows all genes within a genomic window.

```julia
using Pkg
Pkg.add(["GeneticsMakie", "CairoMakie", "DataFrames", "Arrow"])

using GeneticsMakie, CairoMakie, DataFrames, Arrow
url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/gencode.v39lift37.annotation.gtf.gz"
gencode = Arrow.Table("data/gencode/$(splitext(basename(url))[1]).arrow")|> DataFrame
```

We can focus on _CACNA1G_ gene as an example. We can use `GeneticsMakie.plotgenes!` 
to plot all genes within a given chromosome and ± 100 Kb window around gene start and stop sites. 
We can then use `GeneticsMakie.labelgenome` to label the genomic range.

```julia
isdir("figs") || mkdir("figs")
set_theme!(font = "Arial")

gene = "CACNA1G"
chr, start, stop = GeneticsMakie.findgene(gene, gencode)
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
ranges = [start - 1e5, stop + 1e5]
rs = GeneticsMakie.plotgenes!(ax, chr, ranges[1], ranges[2], gencode)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, ranges[1], ranges[2])
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/CACNA1G-gene.png)

We can adjust the height of exons using the `height` keyword argument.
```julia
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs = GeneticsMakie.plotgenes!(ax, chr, ranges[1], ranges[2], gencode; height = 0.1)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, ranges[1], ranges[2])
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/CACNA1G-gene-height0.1.png)

We can change the color of genes and text using the `genecolor` and `textcolor` keyword arguments, respectively.
```julia
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs = GeneticsMakie.plotgenes!(ax, chr, ranges[1], ranges[2], gencode; height = 0.1, genecolor = :mediumorchid3, textcolor = :forestgreen)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, ranges[1], ranges[2])
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/CACNA1G-gene-height0.1-color.png)

Alternatively, we can visualize this locus by passing `gene` as a positional argument and 
`window` as a keyword argument.
```julia
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs = GeneticsMakie.plotgenes!(ax, gene, gencode; window = 1e5, height = 0.1, genecolor = :brown3)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, ranges[1], ranges[2])
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/CACNA1G-gene2.png)

There is no limit to the number of genes we can visualize. Below we visualize a larger 
± 2 Mb window around `gene`.
```julia
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs = GeneticsMakie.plotgenes!(ax, gene, gencode; window = 2e6, height = 0.1)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 2e6, stop + 2e6)
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/CACNA1G-gene3.png)

!!! note "Gene density"
    As some regions have higher gene density than the others, it would be wise (for
    publication purpose) to visualize a smaller genomic window for such gene-dense regions.
    One can instead visualize only the genes of interest (e.g. protein coding genes) by subsetting
    GENCODE accordingly to such features.

We can highlight a gene or sets of genes as below. This can be useful when highlighting genes 
by certain properties such as those that are protein coding or those that are 
[loss-of-function intolerant](https://gnomad.broadinstitute.org/) or those that are 
significant in some sort of gene-level association. 

```julia
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs = GeneticsMakie.plotgenes!(ax, chr, start, (gene, :brown3), gencode; window = 1e5, height = 0.1)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 1e5, start + 1e5)
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/CACNA1G-gene-highlight.png)

```julia
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs = GeneticsMakie.plotgenes!(ax, chr, start, ([gene, "EPN3"], [:brown3, :forestgreen]), gencode; window = 1e5, height = 0.1)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 1e5, start + 1e5)
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/CACNA1G-gene-highlights.png)

Finally, we can make additional modifications on top of the figure as needed using
[__Makie.jl__](https://makie.juliaplots.org/stable/).

```julia
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs = GeneticsMakie.plotgenes!(ax, chr, start, gencode; window = 1e5, height = 0.1)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 1e5, start + 1e5)
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
vlines!(ax, start, color = (:gold, 0.5), linewidth = 0.5)
vlines!(ax, stop, color = (:gold, 0.5), linewidth = 0.5)
f
```
![](../figs/CACNA1G-gene-line.png)

Then we can save the figure as below.
```julia
save("figs/$(gene)-gene.png", f, px_per_unit = 4)
```