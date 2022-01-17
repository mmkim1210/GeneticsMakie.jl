```@meta
CurrentModule = GeneticsMakie
```

# GeneticsMakie

The goal of [__GeneticsMakie__](https://github.com/mmkim1210/GeneticsMakie.jl) is to 
permit data visualization and seamless exploratory data analysis of the human genome within
the larger Julia data science and [__OpenMendel__](https://github.com/OpenMendel) ecosystems.
The package provides convenient wrapper functions for wrangling genetic association results and 
plotting them using [__Makie__](https://makie.juliaplots.org/stable/). Every component of a figure 
can be easily customized and extended. The package generates high-quality, publication-ready figures. 

!!! tip "Getting started"
    Please peruse the documentations of 
    [__Makie__](https://makie.juliaplots.org/stable/), 
    [__CSV__](https://csv.juliadata.org/stable/), 
    [__DataFrames__](https://dataframes.juliadata.org/stable/), 
    and [__SnpArrays__](https://openmendel.github.io/SnpArrays.jl/latest/). 
    Familiarity with these packages allow visualization of any type of genetic and genomic data. 
    [__Makie__](https://makie.juliaplots.org/stable/)'s default layout tools are particularly useful for 
    plotting different genetic and genomic data modalities on separate layers.

!!! note "An usage case"
    If you have run a genome-wide association study (GWAS) at the variant-level, 
    and you would like to eyeball genome-wide significant loci across hundreds of
    phenotypes, then you are in the right place.