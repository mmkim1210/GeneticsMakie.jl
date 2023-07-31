# Munging summary statistics
GWAS summary statistics come in a variety of shapes and flavors, so harmonizing them
is crucial in making our lives easier when trying to visualize their results. We can 
take a peak at GWAS results for height and weight, the two classic anthropometric traits. 

```julia
using Pkg
Pkg.add(["GeneticsMakie", "CSV", "DataFrames", "Arrow"])
```

```julia
using GeneticsMakie, CSV, DataFrames, Arrow, Downloads
gwas = Dict(
    "height" => "https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz",
    "weight" => "https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz"
)
isdir("data/gwas") || mkdir("data/gwas")
dfs = DataFrame[]
for key in keys(gwas)
    url = gwas[key]
    isfile("data/gwas/$(basename(url))") || Downloads.download(url, "data/gwas/$(basename(url))")
    push!(dfs, CSV.read("data/gwas/$(basename(url))", DataFrame; comment = "##", missingstring = ["NA"]))
end
```

To harmonize summary statistics, we run a single command `GeneticsMakie.mungesumstats!`. 
For all downstream plotting functions, we note that GWAS summary statistics should be pre-harmonized.
```julia
GeneticsMakie.mungesumstats!(dfs)
```

To find GWAS loci for each phenotype that are separated by at least 1 Mb, 
we can execute the following command.
```julia
GeneticsMakie.findgwasloci(dfs[1])
GeneticsMakie.findgwasloci(dfs[2])
```

To find every GWAS loci across multiple phenotypes that are separated by at least 1 Mb,
we can instead run the following command.
```julia
GeneticsMakie.findgwasloci(dfs)
```

Such an exhaustive list of loci can then be iterated through and visualized by 
[Plotting LocusZooom](@ref).

We can also find the closest gene to each index SNP in GWAS loci after [Parsing GENCODE](@ref).
```julia
file = "gencode.v39lift37.annotation.gtf.arrow"
gencode = Arrow.Table("data/gencode/$(file)")|> DataFrame

loci = GeneticsMakie.findgwasloci(dfs[1])
GeneticsMakie.findclosestgene(loci, gencode)
GeneticsMakie.findclosestgene(loci, gencode; start = true) # closest gene from gene start site
GeneticsMakie.findclosestgene(loci, gencode; proteincoding = true) # closest "protein-coding" gene
```

To reduce memory intake, we can store and load GWAS summary statistics as Arrow files. 
```julia
for (i, key) in enumerate(keys(gwas))
    Arrow.write("data/gwas/$(key).arrow", dfs[i])
end
```

The summary statistics we are using are on the genome assembly GRCh37. However, newer 
GWAS and annotations may be shared on the latest genome assembly GRCh38. The genomic 
coordinates between the two builds are incompatible; in order to visualize data across 
builds, we need to "liftover" the coordinates all to one genome assembly. GeneticsMakie 
provides functions for perfoming liftover on summary statistics.  

To perform liftover, we need to read a chain file describing how genomic coordinates 
map between one build to the other. Similar to `GeneticsMakie.mungesumstats!`, chromosome 
names are of type `String` and are stripped of the prefix "chr".
```julia
isdir("data/chain") || mkdir("data/chain")
url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
isfile("data/chain/$(basename(url))") || Downloads.download(url, "data/chain/$(basename(url))")

chain = GeneticsMakie.readchain("data/chain/$(basename(url))")
```

Additionally, we will need to read in FASTA files describing both the reference and 
the query genome builds (GRCh37 and GRCh38 respectively). These are used to properly 
liftover indels that were on the opposite strand. You will need to uncompress and 
index these files to use them with the `FASTX` julia package.

```julia
isdir("data/fasta") || mkdir("data/fasta")
urls = [
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
]
for url in urls
    isfile("data/fasta/$(basename(url))") || Downloads.download(url, "data/fasta/$(basename(url))")
end
```

```bash
gunzip data/fasta/hg19.fa.gz
gunzip data/fasta/hg38.fa.gz
samtools faidx data/fasta/hg19.fa
samtools faidx data/fasta/hg38.fa
```

```julia
import FASTX: FASTA
hg19io = open("data/fasta/hg19.fa")
hg19fa = FASTA.Reader(hg19io)
hg19fai = FASTA.Index("data/fasta/hg19.fa.fai")
FASTA.index!(hg19fa, hg19fai)
hg38io = open("data/fasta/hg38.fa")
hg38fa = FASTA.Reader(hg38io)
hg38fai = FASTA.Index("data/fasta/hg38.fa.fai")
FASTA.index!(hg38fa, hg38fai)
```

With the chain file and the FASTA files loaded, we can now perform liftover on our 
GWAS. `GeneticsMakie.liftoversumstats!` will liftover the sumstats in place and return 
a NamedTuple of `unmapped`) unmapped variants still on the original build and `multiple`) 
variants on the target build that has mapped to multiple positions. Liftover will 
take a while to complete, especially summary statistics with many variants.
```julia
dfs_hg38 = deepcopy(dfs)
unmapped, multiple = GeneticsMakie.liftoversumstats!(dfs_hg38, chain, hg19fa, hg38fa;
                                                     multiplematches = :warning,
                                                     pickreference = :longest,
                                                     indelref = :start,
                                                     extendambiguous = true)
close(hg19io)
close(hg38io)
```

