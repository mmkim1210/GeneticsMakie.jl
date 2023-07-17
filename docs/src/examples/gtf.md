# Parsing GENCODE

Install the relevant packages in the usual way.
```julia
using Pkg
Pkg.add(["GeneticsMakie", "CSV", "DataFrames", "Arrow"])
```

To plot genes and isoforms, we need a transcriptome annotation. We can use 
the latest [GENCODE](https://www.gencodegenes.org/human/) annotation for 
the human genome (GRCh37), where we download the comprehensive 
gene annotation file in [GTF format](https://www.gencodegenes.org/pages/data_format.html).
We recommend having at least 16 GB RAM available for loading GENCODE annotation.

```julia
using GeneticsMakie, CSV, DataFrames, Arrow, Downloads
isdir("data") || mkdir("data")
url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/gencode.v39lift37.annotation.gtf.gz"
file = basename(url)
isdir("data/gencode") || mkdir("data/gencode")
isfile("data/gencode/$(file)") || Downloads.download(url, "data/gencode/$(file)")
h = ["seqnames", "source", "feature", "start", "end", "score", "strand", "phase", "info"]
gencode = CSV.read("data/gencode/$(file)", DataFrame; delim = "\t", comment = "#", header = h)
```

!!! warning "Human genome build"
    The latest human genome assembly is [GRCh38](https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38.p14), but we use an annotation with coordinates 
    from the older version (GRCh37), because a lot of the GWAS results are shared in 
    GRCh37 genomic coordinates. Make sure to either use the matching human genome 
    build or perform liftover on the GWAS summary statistics (as described in 
    [Munging summary statistics](@ref)) when visualizing your results. 

The ninth column of a [GTF file](https://uswest.ensembl.org/info/website/upload/gff.html) 
contains rich information about features, so we can parse this column.
```julia
GeneticsMakie.parsegtf!(gencode)
```

!!! info "Chromosome names"
    Chromosome names are munged to not contain “chr” prefix, and their type is `String`,
    since there could be non-numerical chromosome names, such as sex chromosomes and mitochondrial genome.

To reduce memory intake, we can also subset `gencode` to most commonly used columns
in downstream analyses.
```julia
select!(gencode, :seqnames, :feature, :start, :end, :strand, :gene_id, :gene_name, :gene_type, :transcript_id)
```

To further reduce memory intake, we can instead store and load GENCODE annotation as an Arrow file. 
```julia
Arrow.write("data/gencode/$(splitext(file)[1]).arrow", gencode)
gencode = Arrow.Table("data/gencode/$(splitext(file)[1]).arrow")|> DataFrame
```

Other transcriptome annotations, such as one from RefSeq, can be used for plotting functions 
as long as they contain the above columns with the right column names.

Once `gencode` is ready, we can look up where a gene is on the human genome.
```julia
GeneticsMakie.findgene("RBFOX1", gencode)
GeneticsMakie.findgene("ENSG00000078328", gencode)
```

!!! tip "Gene names"
    Make sure to use the correct gene name in case the gene cannot be found.
    The latest gene names can be looked up in databases such as [GeneCards](https://www.genecards.org/).
