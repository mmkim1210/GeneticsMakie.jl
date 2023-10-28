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
include("mungesumstats.jl")
include("plotlocus.jl")
include("plotloops.jl")
include("plotqq.jl")
include("plotgwas.jl")
include("gwas.jl")
include("genome.jl")
include("plotld.jl")
include("plotrg.jl")

end