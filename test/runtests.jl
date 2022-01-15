using GeneticsMakie
using Test
using CairoMakie
using Makie.GeometryBasics
using CSV
using DataFrames
using SnpArrays
using Statistics

@testset "Plotting genes/isoforms" begin
    gencode = CSV.read("data/gencode.gtf", DataFrame)
    GeneticsMakie.parsegtf!(gencode)
    @test ncol(gencode) == 14

    gene = "KMT2E"
    chr, start, stop = GeneticsMakie.findgene(gene, gencode)
    @test ("7", 104581390, 104755466) == (chr, start, stop)
    range1 = start - 1e6
    range2 = stop + 1e6

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, range1, range2, gencode)
    GeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, range1, range2)
    save("gene.png", f)
    @test isfile("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, start, gencode)
    save("gene.png", f)
    @test isfile("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, gene, gencode)
    save("gene.png", f)
    @test isfile("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, range1, range2, (["KMT2E", "SETD1A"], ["red", "blue"]), gencode)
    save("gene.png", f)
    @test isfile("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, start, (["KMT2E", "SETD1A"], ["red", "blue"]), gencode)
    save("gene.png", f)
    @test isfile("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, gene, (["KMT2E", "SETD1A"], ["red", "blue"]), gencode)
    save("gene.png", f)
    @test isfile("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, range1, range2, ("KMT2E", "red"), gencode)
    save("gene.png", f)
    @test isfile("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, start, ("KMT2E", "red"), gencode)
    save("gene.png", f)
    @test isfile("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, gene, ("KMT2E", "red"), gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    # plotisoforms
end

@testset "Munging summmary stats" begin
end

@testset "Plotting LD" begin
    LD = rand(10, 10)
    LD = LD + transpose(LD)
    for i in 1:size(LD, 1)
        LD[i, i] = 1
    end
    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotld!(ax, LD)
    save("ld.png", f)
    @test isfile("ld.png")
    rm("ld.png")
end

@testset "Plotting LocusZoom" begin
end

@testset "Plotting GWAS" begin
end

@testset "Plotting QQ plot" begin
end

@testset "Plotting TWAS" begin
end

@testset "Plotting correlation" begin
end