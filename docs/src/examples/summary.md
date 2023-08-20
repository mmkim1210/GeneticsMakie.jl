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

Sometimes, it will be useful or necessary to harmonize the variant names as well. 
One set of summary statistics may have the variants named by their chromosome and 
position while another set may use the rsIDs from dbSNP to label the SNPs. With the 
following code, we can harmonize variant names for our summary statistics to a specific 
dbSNP version, falling back to `"$CHR:$BP:$A1:$A2"` if the variant is absent from 
the dbSNP database.  

First, we will need to normalize our summary statistics to the reference genome. Sometimes 
summary statistics may have the alleles swapped by certain tools, so we will normalize 
the summary statistics to match the GRCh37 reference build, which the dbSNP database 
is also matched to.  

The reference build is stored in the FASTA format. We will need version GRCh37 as 
that is what our summary statistics are based on. UCSC shares these reference genomes 
freely on their site; we will download the files from there and to uncompress and 
index these files to use them with the `FASTX` julia package. We also download GRCh38 
here too, as it will be relevant for our later examples (but not for this particular 
example of variant name harmonization).

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
hg19 = FASTA.Reader(open("data/fasta/hg19.fa"), index = "data/fasta/hg19.fa.fai")
GeneticsMakie.normalizesumstats!(dfs, hg19)
close(hg19)
```

Now, we will harmonize the variant names to dbSNP 151. Several versions of this reference 
can be downloaded from the NIH site; we will download the version containing common 
variants.
```julia
isdir("data/vcf") || mkdir("data/vcf")
urls = [
    "https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz",
]
for url in urls
    isfile("data/vcf/$(basename(url))") || Downloads.download(url, "data/vcf/$(basename(url))")
end
```

There are two different ways to use the `harmonizevariantnames!` function; either 
you load in the VCF file into memory as a DataFrame or you use the CSV.Rows iterator 
to read in one row at a time. Loading the VCF file into memory will speed up repeated 
calls of `harmonizevariantnames!` at the cost of high memory usage while iterating 
row by row will keep memory usage to a minimum but necessitate iterating through the 
whole VCF file every function call. Additionally, when iterating, the file is expected 
to be sorted by chromosome and position (i.e. `bcftools sort`); the dbSNP VCF files 
are already sorted.

```julia
dbsnp_dataframe = CSV.read("data/vcf/00-common_all.vcf.gz", DataFrame; delim = "\t", comment = "##")
GeneticsMakie.harmonizevariantnames!(dfs, dbsnp_dataframe)

dbsnp_iterator = CSV.Rows("data/vcf/00-common_all.vcf.gz", delim = "\t", comment = "##")
GeneticsMakie.harmonizevariantnames!(dfs, dbsnp_iterator)
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
hg19 = FASTA.Reader(open("data/fasta/hg19.fa"), index = "data/fasta/hg19.fa.fai")
hg38 = FASTA.Reader(open("data/fasta/hg38.fa"), index = "data/fasta/hg38.fa.fai")
```

With the chain file and the FASTA files loaded, we can now perform liftover on our 
GWAS. `GeneticsMakie.liftoversumstats!` will liftover the sumstats in place and return 
a DataFrame of `unmapped`) unmapped variants still on the original build and `multiple`) 
variants on the target build that has mapped to multiple positions. Liftover will 
take a while to complete, especially summary statistics with many variants.
```julia
dfs_hg38 = deepcopy(dfs)
unmapped, multiple = GeneticsMakie.liftoversumstats!(dfs_hg38, hg19fa, hg38fa, chain; 
                                                     multiplematches = :warning,
                                                     whichreference = :first_silent,
                                                     indelreference = :start,
                                                     extendambiguous = true)
close(hg19)
close(hg38)
```

