module GeneticsMakie

using CairoMakie
using Makie.GeometryBasics
using DataFrames
using SnpArrays
using Statistics
using Distributions

include("parsegtf.jl")
include("plotgenes.jl")
include("plotisoforms.jl")
include("plotld.jl")
include("mungesumstats.jl")
include("plotlocus.jl")
include("plotqq.jl")
include("plotgwas.jl")
include("plotrg.jl")
include("gwas.jl")
include("genome.jl")

end