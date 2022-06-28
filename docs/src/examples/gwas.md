# Plotting GWAS

After [Munging summary statistics](@ref), we can use `GeneticsMakie.plotgwas!` 
to draw Manhattan plots.

```julia
using Pkg
Pkg.add(["GeneticsMakie", "CairoMakie", "DataFrames", "Arrow"])
```

```julia
using GeneticsMakie, CairoMakie, DataFrames, Arrow
dfs = DataFrame[]
for key in ["height", "weight"]
    push!(dfs, Arrow.Table("data/gwas/$(key).arrow")|> DataFrame)
end
titles = ["Height (Yengo et al. 2018)", "Weight (Yengo et al. 2018)"]
```

```julia
f = Figure(resolution = (408, 792))
axs = [Axis(f[i, 1]) for i in 1:length(titles)]
for i in eachindex(titles)
    GeneticsMakie.plotgwas!(axs[i], dfs[i])
    hidespines!(axs[i], :t, :r)
    Label(f[i, 1, Top()], text = "$(titles[i])", textsize = 8)
    rowsize!(f.layout, i, 50)
    i == length(titles) ? axs[i].xlabel = "Chromosome" : axs[i].xlabel = ""
end
rowgap!(f.layout, 10)
resize_to_layout!(f)
f
```
![](../figs/manhattan.png)

By default, `GeneticsMakie.plotgwas!` highlights the genome-wide significant threshold
and corresponding significant variants. We can turn off this option by using the `linecolor` and `scattercolor` keyword arguments. 

```julia
f = Figure(resolution = (408, 792))
axs = [Axis(f[i, 1]) for i in 1:length(titles)]
for i in eachindex(titles)
    GeneticsMakie.plotgwas!(axs[i], dfs[i]; linecolor = nothing, scattercolor = nothing)
    hidespines!(axs[i], :t, :r)
    Label(f[i, 1, Top()], text = "$(titles[i])", textsize = 8)
    rowsize!(f.layout, i, 50)
    i == length(titles) ? axs[i].xlabel = "Chromosome" : axs[i].xlabel = ""
end
rowgap!(f.layout, 10)
resize_to_layout!(f)
f
```
![](../figs/manhattan-nocolor.png)

We can color even and odd chromosomes with different colors by using the `chromcolors` keyword argument. 

```julia
f = Figure(resolution = (408, 792))
axs = [Axis(f[i, 1]) for i in 1:length(titles)]
for i in eachindex(titles)
    GeneticsMakie.plotgwas!(axs[i], dfs[i]; linecolor = nothing, scattercolor = nothing, 
        chromcolors = ["#389826", "#9658B2"])
    hidespines!(axs[i], :t, :r)
    Label(f[i, 1, Top()], text = "$(titles[i])", textsize = 8)
    rowsize!(f.layout, i, 50)
    i == length(titles) ? axs[i].xlabel = "Chromosome" : axs[i].xlabel = ""
end
rowgap!(f.layout, 10)
resize_to_layout!(f)
f
```
![](../figs/manhattan-chromosomes.png)