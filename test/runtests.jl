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
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, "6", range1, range2, gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, start, gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, gene, gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, range1, range2, (["KMT2E", "SETD1A"], ["red", "blue"]), gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, start, (["KMT2E", "SETD1A"], ["red", "blue"]), gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, gene, (["KMT2E", "SETD1A"], ["red", "blue"]), gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, range1, range2, ("KMT2E", "red"), gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, "6", range1, range2, ("KMT2E", "red"), gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, chr, start, ("KMT2E", "red"), gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgenes!(ax, gene, ("KMT2E", "red"), gencode)
    save("gene.png", f)
    @test isfile("gene.png")
    rm("gene.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotisoforms!(ax, gene, gencode; orderby = ["ENST00000667857", "ENST00000482560", "Random"])
    save("isoform.png", f)
    @test isfile("isoform.png")
    rm("isoform.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotisoforms!(ax, gene, gencode; text = :l)
    save("isoform.png", f)
    @test isfile("isoform.png")
    rm("isoform.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotisoforms!(ax, gene, gencode; text = :b)
    save("isoform.png", f)
    @test isfile("isoform.png")
    rm("isoform.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotisoforms!(ax, gene, gencode; text = :r)
    save("isoform.png", f)
    @test isfile("isoform.png")
    rm("isoform.png")
end

@testset "Munging summmary stats" begin

end

@testset "Plotting GWAS / QQ" begin
    gwas = DataFrame(CHR = rand(string.(collect(1:22)), 1000),
        BP = rand(1:1000, 1000),
        P = rand(1000))

    f = Figure()
    ax = Axis(f[1, 1])
    coord, ymaxs, xmax, ticks = GeneticsMakie.coordinategwas([gwas])
    GeneticsMakie.plotgwas!(ax, coord, 1, ymaxs[1], xmax, ticks)
    save("manhattan.png", f)
    @test isfile("manhattan.png")
    rm("manhattan.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotqq!(ax, gwas)
    save("qq.png", f)
    @test isfile("qq.png")
    rm("qq.png")
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
    kgp = SnpData("data/kgp")
    gwas = CSV.read("data/locus.csv", DataFrame)

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotlocus!(ax, "7", 104581390, 104755466, gwas)
    save("locuszoom.png", f)
    @test isfile("locuszoom.png")
    rm("locuszoom.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotlocus!(ax, "7", 104755466, 104581390, gwas)
    save("locuszoom.png", f)
    @test isfile("locuszoom.png")
    rm("locuszoom.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotlocus!(ax, "7", 104581390, 104755466, gwas; ld = kgp)
    save("locuszoom.png", f)
    @test isfile("locuszoom.png")
    rm("locuszoom.png")    

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotlocus!(ax, "7", 104581390, 104755466, gwas; ld = (kgp, "rs11764361"))
    save("locuszoom.png", f)
    @test isfile("locuszoom.png")
    rm("locuszoom.png")    

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotlocus!(ax, "7", 104581390, 104755466, gwas; ld = (kgp, "rs111931861"))
    save("locuszoom.png", f)
    @test isfile("locuszoom.png")
    rm("locuszoom.png")
end

@testset "Plotting QQ plot" begin
    P = rand(1000)
    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotqq!(ax, P)
    save("qq.png", f)
    @test isfile("qq.png")
    rm("qq.png")

    df = DataFrame(P = P)
    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotqq!(ax, P)
    save("qq.png", f)
    @test isfile("qq.png")
    rm("qq.png")
end

@testset "Plotting correlation" begin
    f = Figure()
    ax = Axis(f[1, 1])
    n = 10
    GeneticsMakie.plotrg!(ax, (rand(n, n) .* 2 .- 1), string.(1:n), circle = true)
    colsize!(f.layout, 1, Aspect(1, 1))
    rowsize!(f.layout, 1, 18 * n)
    save("cor.png", f)
    @test isfile("cor.png")
    rm("cor.png")
end

@testset "Plotting TWAS" begin
end