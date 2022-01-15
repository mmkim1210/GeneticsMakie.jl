using GeneticsMakie
using Test
using CSV
using CairoMakie
using DataFrames
using SnpArrays
using Statistics

@testset "Parsing GENCODE" begin
    gencode = CSV.read("data/gencode.gtf", DataFrame)
    GeneticsMakie.parsegtf!(gencode)
    @test ncol(gencode) == 14
end

@testset "Plotting genes" begin
end

@testset "Plotting isoforms" begin
end

@testset "Munging summmary stats" begin
end

@testset "Plotting LD" begin
end

@testset "Plotting LocusZoom" begin
end

@testset "Plotting GWAS" begin
    # @test isfile("gwas.png")
    # rm("gwas.png")
end

@testset "Plotting qq" begin
end

@testset "Plotting TWAS" begin
end