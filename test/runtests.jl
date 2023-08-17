using GeneticsMakie
using Test
using CairoMakie
using CSV
using DataFrames
import FASTX: FASTA
using SnpArrays
using Statistics

@testset "Plotting genes/isoforms" begin
    gencode = CSV.read("data/gencode.gtf", DataFrame)
    GeneticsMakie.parsegtf!(gencode)
    @test ncol(gencode) == 14

    gene = "KMT2E"
    chr, start, stop = GeneticsMakie.findgene(gene, gencode)
    @test ("7", 104581390, 104755466) == (chr, start, stop)
    _, start, stop = GeneticsMakie.findgene("ENSG00000005483", gencode)
    @test (104581390, 104755466) == (start, stop)
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
    gwas = CSV.read("data/sumstats.csv", DataFrame)
    GeneticsMakie.mungesumstats!(gwas)
    @test ncol(gwas) == 7

    gwas = CSV.read("data/sumstats.csv", DataFrame)
    select!(gwas, Not("BETA"))
    GeneticsMakie.mungesumstats!(gwas)
    @test ncol(gwas) == 6

    gwas = CSV.read("data/sumstats.csv", DataFrame)
    select!(gwas, "#CHROM", :POS, :PVAL, :BETA => (col -> exp.(col)) => :OR, :SE)
    GeneticsMakie.mungesumstats!(gwas)
    @test ncol(gwas) == 5

    gwas = CSV.read("data/sumstats.csv", DataFrame)
    select!(gwas, "#CHROM", :POS, :PVAL, :BETA)
    GeneticsMakie.mungesumstats!(gwas)
    @test ncol(gwas) == 5

    gwas = CSV.read("data/sumstats.csv", DataFrame)
    select!(gwas, "#CHROM", :POS, :PVAL, :BETA => (col -> exp.(col)) => :OR)
    GeneticsMakie.mungesumstats!(gwas)
    @test ncol(gwas) == 5

    gwas = CSV.read("data/sumstats.csv", DataFrame)
    select!(gwas, "#CHROM", :POS, :PVAL)
    GeneticsMakie.mungesumstats!(gwas)
    @test ncol(gwas) == 4

    df = GeneticsMakie.findgwasloci(gwas)
    @test nrow(df) == 1
    df = GeneticsMakie.findgwasloci([gwas, gwas])
    @test nrow(df) == 1

    chain = GeneticsMakie.readchain("data/hg19ToHg38.over.chain.gz")
    @test nrow(chain) == 53950
    @test ncol(chain) == 11

    run(`gzip -d data/GRCh37.p13.chr20.fa.gz`)
    run(`gzip -d data/GRCh38.p13.chr20.fa.gz`)
    fasta37 = FASTA.Reader(open("data/GRCh37.p13.chr20.fa"), index = "data/GRCh37.p13.chr20.fa.fai")
    fasta38 = FASTA.Reader(open("data/GRCh38.p13.chr20.fa"), index = "data/GRCh38.p13.chr20.fa.fai")

    gwas = CSV.read("data/sumstats.w_indels.csv.gz", DataFrame)
    GeneticsMakie.mungesumstats!(gwas)
    unmapped, multiple =
    GeneticsMakie.liftoversumstats!([gwas], fasta37, fasta38, chain; referenceorder = false)
    @test nrow(unmapped[1]) == 0
    @test nrow(multiple[1]) == 0

    gwas = CSV.read("data/sumstats.w_indels.csv", DataFrame)
    GeneticsMakie.mungesumstats!(gwas)
    unmapped, multiple =
    GeneticsMakie.liftoversumstats!([gwas], fasta37, fasta38, chain; referenceorder = true)
    bcftools_score = CSV.read("data/adhd_Demontis.38.vcf", DataFrame, comment = "##")
    select!(bcftools_score,
            :ID => :SNP, "#CHROM" => (col -> replace.(col, "chr" => "")) => :CHR,
            :POS => :BP, :REF => :A1, :ALT => :A2)
    @test nrow(innerjoin(gwas, bcftools_score, on = [:SNP, :CHR, :BP, :A1, :A2])) == 1479
    @test nrow(unmapped[1]) == 0
    @test nrow(multiple[1]) == 0

    close(fasta37)
    close(fasta38)
    run(`gzip data/GRCh37.p13.chr20.fa`)
    run(`gzip data/GRCh38.p13.chr20.fa`)

    push!(chain, chain[nrow(chain), :])
    @test_throws ErrorException GeneticsMakie.findnewcoord("Y",
                                                           59282592,
                                                           chain;
                                                           multiplematches = :error)
    @test length(GeneticsMakie.findnewcoord("Y", 59282592, chain;
                                            multiplematches = :warning)) == 2
end

@testset "Plotting GWAS" begin
    gwas = DataFrame(
        CHR = rand(string.(collect(1:22)), 1000),
        BP = rand(1:1000, 1000),
        P = rand(1000)
    )

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotgwas!(ax, gwas)
    save("manhattan.png", f)
    @test isfile("manhattan.png")
    rm("manhattan.png")
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
    gwas.CHR = string.(gwas.CHR)

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotlocus!(ax, "7", 103581390, 105755466, gwas)
    save("locuszoom.png", f)
    @test isfile("locuszoom.png")
    rm("locuszoom.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotlocus!(ax, "7", 105755466, 103581390, gwas)
    save("locuszoom.png", f)
    @test isfile("locuszoom.png")
    rm("locuszoom.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotlocus!(ax, "7", 103581390, 105755466, gwas; ld = kgp)
    save("locuszoom.png", f)
    @test isfile("locuszoom.png")
    rm("locuszoom.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotlocus!(ax, "7", 103581390, 105755466, gwas; ld = (kgp, "rs11764361"))
    save("locuszoom.png", f)
    @test isfile("locuszoom.png")
    rm("locuszoom.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotlocus!(ax, "7", 103581390, 105755466, gwas; ld = (kgp, "rs111931861"))
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

@testset "Plotting loops" begin
    loopdf = CSV.read("data/loops.csv", DataFrame)
    loopdf.chr1 = string.(loopdf.chr1)
    loopdf.chr2 = string.(loopdf.chr2)

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotloops!(ax, "7", 103581390, 105755466, loopdf)
    save("loops.png", f)
    @test isfile("loops.png")
    rm("loops.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotloops!(ax, "7", 105755466, 103581390, loopdf)
    save("loops.png", f)
    @test isfile("loops.png")
    rm("loops.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotloops!(ax, "7", 103581390, 105755466, loopdf; ymax = 150)
    save("loops.png", f)
    @test isfile("loops.png")
    rm("loops.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotloops!(ax, "7", 103581390, 105755466, loopdf;
                             linewidth = 10)
    save("loops.png", f)
    @test isfile("loops.png")
    rm("loops.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotloops!(ax, "7", 103581390, 105755466, loopdf;
                             colorarc = :royalblue)
    save("loops.png", f)
    @test isfile("loops.png")
    rm("loops.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotloops!(ax, "7", 103581390, 105755466, loopdf;
                             colorend = :royalblue)
    save("loops.png", f)
    @test isfile("loops.png")
    rm("loops.png")

    f = Figure()
    ax = Axis(f[1, 1])
    GeneticsMakie.plotloops!(ax, "7", 103581390, 105755466, loopdf;
                             resolution = 10000)
    save("loops.png", f)
    @test isfile("loops.png")
    rm("loops.png")
end
