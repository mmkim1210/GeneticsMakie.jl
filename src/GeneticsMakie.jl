module GeneticsMakie

using Makie
using CairoMakie
using Makie.GeometryBasics
using DataFrames
using Chain
using DataFramesMeta
using CSV
using Colors
using ColorSchemes
using SnpArrays
using Statistics
using Distributions

include("parsegtf.jl")
include("plotgenes.jl")
include("plotisoforms.jl")
include("plotld.jl")
include("plotlocus.jl")
# include("mungesumstats.jl")
# include("plotqq.jl")
# include("plotgwas.jl")
# include("plotpca.jl")
# include("plotadmixture.jl")
# include("plotwas.jl")
# include("plotrg.jl")
# include("plotcoloc.jl")

end