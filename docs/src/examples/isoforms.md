# Plotting isoforms

After [Parsing GENCODE](@ref), it is possible to plot isoform bodies. We can focus on _NRXN1_ gene as our initial example. `GeneticsMakie.plotisoforms!` returns 
genomic coordinates for the gene of interest so that an appropriate label can be passed onto 
`GeneticsMakie.labelgenome`. _NRXN1_ gene has many isoforms as we see below, and even more 
isoforms are likely to be discovered in the future. For this reason, plotting isoforms of multiple genes is not available. 

```julia
using Pkg
Pkg.add(["GeneticsMakie", "CairoMakie", "DataFrames", "Arrow"])

using GeneticsMakie, CairoMakie, DataFrames, Arrow
url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/gencode.v39lift37.annotation.gtf.gz"
gencode = Arrow.Table("data/gencode/$(splitext(basename(url))[1]).arrow")|> DataFrame
```

```julia
isdir("figs") || mkdir("figs")
set_theme!(font = "Arial")

gene = "NRXN1"
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs, chr, range1, range2 = GeneticsMakie.plotisoforms!(ax, gene, gencode; height = 0.1)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, range1, range2)
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/NRXN1-isoform.png)

To save some space, we can plot the isoform labels on the left by using the `text` keyword argument. 

```julia
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs, chr, range1, range2 = GeneticsMakie.plotisoforms!(ax, gene, gencode; height = 0.1, text = :l)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, range1, range2)
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/NRXN1-isoform-left.png)

We can change the color of isoforms and text using the `isoformcolor` and `textcolor` keyword arguments, respectively.
We look at a different gene, _GRIN2A_.

```julia
gene = "GRIN2A"
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs, chr, range1, range2 = GeneticsMakie.plotisoforms!(ax, gene, gencode; isoformcolor = :forestgreen, height = 0.1)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, range1, range2)
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/GRIN2A-isoform-color.png)

We can change the order of isoforms by using the `orderby` keyword argument.

```julia
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs, chr, range1, range2 = GeneticsMakie.plotisoforms!(ax, gene, gencode; orderby = ["ENST00000675189", "ENST00000463531"], isoformcolor = :forestgreen, height = 0.1)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, range1, range2)
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/GRIN2A-isoform-order.png)

```julia
f = Figure(size = (306, 792))
ax = Axis(f[1, 1])
rs, chr, range1, range2 = GeneticsMakie.plotisoforms!(ax, gene, gencode; orderby = ["ENST00000675189", "ENST00000463531"], 
    highlight = (["ENST00000675398", "ENST00000562109"], fill(:brown3, 2)), isoformcolor = :forestgreen, height = 0.1)
GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, range1, range2)
rowsize!(f.layout, 1, rs)
resize_to_layout!(f)
f
```
![](../figs/GRIN2A-isoform-highlights.png)