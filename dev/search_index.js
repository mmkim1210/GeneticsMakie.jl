var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [GeneticsMakie]","category":"page"},{"location":"api/#GeneticsMakie.coordinategwas-Tuple{Vector{DataFrames.DataFrame}}","page":"API","title":"GeneticsMakie.coordinategwas","text":"coordinategwas(gwas::Vector{DataFrame}; freey::Bool)\n\nDetermine shared coordinates for a set of gwas.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.findclosestgene-Tuple{AbstractString, Real, DataFrames.DataFrame}","page":"API","title":"GeneticsMakie.findclosestgene","text":"findclosestgene(chr::AbstractString, bp::Real, gencode::DataFrame; start::Bool, proteincoding::Bool)\nfindclosestgene(df::DataFrame, gencode::DataFrame; start::Bool, proteincoding::Bool)\n\nFind the closest gene(s) to a genomic coordinate or a list of genomic coordinates using gencode. \n\nOptionally, the closest gene can be defined from the gene start site using start, and only protein coding genes can be considered using proteincoding.  The default start and proteincoding are false.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.findgene-Tuple{AbstractString, DataFrames.DataFrame}","page":"API","title":"GeneticsMakie.findgene","text":"findgene(gene::AbstractString, gencode::DataFrame)\n\nFind chromosome, gene start, and gene stop sites for the gene of interest.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.findgwasloci-Tuple{DataFrames.DataFrame}","page":"API","title":"GeneticsMakie.findgwasloci","text":"findgwasloci(gwas::DataFrame; p::Real)\nfindgwasloci(gwas::Vector{DataFrame}; p::Real)\n\nFind genome-wide significant loci for gwas that are separated from each other by at least 1 Mb.\n\nAlternatively, find genome-wide significant loci across multiple gwas that  are all separated by at least 1 Mb. p determines the genome-wide significance threshold,  which is 5e-8 by default.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.labelgenome-Tuple{GridLayoutBase.GridPosition, AbstractString, Real, Real}","page":"API","title":"GeneticsMakie.labelgenome","text":"labelgenome(g::GridPosition, chromosome::AbstractString, range1::Real, range2::Real)\n\nLabel g with a given chromosome and genomic range between range1 and range2.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.mungesumstats!-Tuple{Vector{DataFrames.DataFrame}}","page":"API","title":"GeneticsMakie.mungesumstats!","text":"mungesumstats!(gwas::DataFrame)\nmungesumstats!(gwas::Vector{DataFrame})\n\nMunge gwas by harmonizing the names of columns, their types, and P values, among others.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.parsegtf!-Tuple{DataFrames.DataFrame}","page":"API","title":"GeneticsMakie.parsegtf!","text":"parsegtf!(gencode::DataFrame)\n\nParse gencode by extracting gene_id, gene_name, gene_type, transcript_id, transcript_support_level information from the info column.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.plotgenes!-Tuple{Makie.MakieLayout.Axis, AbstractString, Real, Real, DataFrames.DataFrame}","page":"API","title":"GeneticsMakie.plotgenes!","text":"plotgenes!(ax::Axis, chromosome::AbstractString, range1::Real, range2::Real, gencode::DataFrame; kwargs)\nplotgenes!(ax::Axis, chromosome::AbstractString, bp::Real, gencode::DataFrame; kwargs)\nplotgenes!(ax::Axis, gene::AbstractString, gencode::DataFrame; kwargs)\n\nPlot collapsed gene bodies for genes within a given chromosome and genomic range  between range1 and range2.\n\nAlternatively, plot within a given chromosome and a certain window around a  genomic coordinate bp or plot within a certain window around gene.\n\nArguments\n\nheight::Real = 0.25: the height of exons.\ngenecolor = :royalblue: the color of genes.\ntextcolor = :black: the color of gene labels.\nwindow::Real = 1e6: the window around bp or gene.                                       \n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.plotgenes!-Tuple{Makie.MakieLayout.Axis, AbstractString, Real, Real, Tuple{AbstractVector, AbstractVector}, DataFrames.DataFrame}","page":"API","title":"GeneticsMakie.plotgenes!","text":"plotgenes!(ax::Axis, chromosome::AbstractString, range1::Real, range2::Real, highlight::Tuple{AbstractVector, AbstractVector}, gencode::DataFrame; height::Real)\nplotgenes!(ax::Axis, chromosome::AbstractString, bp::Real, highlight::Tuple{AbstractVector, AbstractVector}, gencode::DataFrame; window::Real, height::Real)\nplotgenes!(ax::Axis, gene::AbstractString, highlight::Tuple{AbstractVector, AbstractVector}, gencode::DataFrame; window::Real, height::Real)\n\nPlot gene bodies with a vector of genes highlighted by a vector of colors via highlight.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.plotgwas!-Tuple{Makie.MakieLayout.Axis, DataFrames.DataFrame, Int64, Real, Real, DataFrames.DataFrame}","page":"API","title":"GeneticsMakie.plotgwas!","text":"plotgwas!(ax::Axis, df::DataFrame, i::Int, ymax::Real, xmax::Real, ticks::DataFrame; ystep::Real, sigline::Bool, sigcolor::Bool)\n\nPlot gwas results with coordinates from coordinategwas, namely df for phenotype i with x and y limits, xmax and ymax, respectively. The position of x ticks are determined by ticks.\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.plotisoforms!-Tuple{Makie.MakieLayout.Axis, AbstractString, DataFrames.DataFrame}","page":"API","title":"GeneticsMakie.plotisoforms!","text":"plotisoforms!(ax::Axis, gene::AbstractString, gencode::DataFrame; kwargs)\n\nPlot each isoform of a given gene on a separate row.\n\nArguments\n\norderby::Union{Nothing, AbstractVector{<:AbstractString}} = nothing: the order of isoforms.\nheight::Real = 0.25: the height of exons.\nisoformcolor = :royalblue: the color of isoforms.\ntextcolor = :black: the color of isoform labels.\ntext::Union{Bool, Symbol} = :top: the position of isoform labels. \n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.plotld!-Tuple{Makie.MakieLayout.Axis, AbstractMatrix}","page":"API","title":"GeneticsMakie.plotld!","text":"plotld!(ax::Axis, LD::AbstractMatrix; color::AbstractString)\n\nVisualize a symmetric correlation matrix LD with the diagonal elements on the x axis. color can be either \"red\", \"green\", \"blue\", or \"black\".\n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.plotlocus!-Tuple{Makie.MakieLayout.Axis, AbstractString, Real, Real, DataFrames.DataFrame}","page":"API","title":"GeneticsMakie.plotlocus!","text":"plotlocus!(ax::Axis, chromosome::AbstractString, range1::Real, range2::Real, gwas::DataFrame; kwargs)\nplotlocus!(ax::Axis, chromosome::AbstractString, bp::Real, gwas::DataFrame; kwargs)\nplotlocus!(ax::Axis, gene::AbstractString, gwas::DataFrame, gencode::DataFrame; kwargs)\n\nPlot gwas results within a given chromosome and genomic range between range1  and range2.\n\nAlternatively, plot within a given chromosome and a certain window around a  genomic coordinate bp or plot within a certain window around gene.\n\nArguments\n\nld::Union{Nothing, SnpData, Tuple{SnpData, Union{AbstractString, Tuple{AbstractString, Int}}}} = nothing:    the reference panel for which LD is calculated.\nymax::Real: the maximum value for y axis. \nwindow::Real = 1e6: the window around bp or gene. \n\n\n\n\n\n","category":"method"},{"location":"api/#GeneticsMakie.plotqq!-Tuple{Makie.MakieLayout.Axis, DataFrames.DataFrame}","page":"API","title":"GeneticsMakie.plotqq!","text":"plotqq!(ax::Axis, df::DataFrame; kwargs)\nplotqq!(ax::Axis, P::AbstractVector; kwargs)\n\nPlot QQ plot of P values where the expected distribution is the uniform distribution.\n\nKeyword arguments include xstep::Real and ystep::Real for x and y axes ticks step sizes.\n\n\n\n\n\n","category":"method"},{"location":"examples/genes/#Plotting-genes","page":"Plotting genes","title":"Plotting genes","text":"","category":"section"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"After Parsing GENCODE, we can start plotting gene bodies.  GeneticsMakie.jl is transparent  in that it shows all genes within a genomic window.","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"using Pkg\nPkg.add(\"CairoMakie\")","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"We can focus on CACNA1G gene as an example. We can use GeneticsMakie.plotgenes!  to plot all genes within a given chromosome and ± 100 Kb window around gene start and stop sites.  We can then use GeneticsMakie.labelgenome to label the genomic range.","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"using CairoMakie\nisdir(\"figs\") || mkdir(\"figs\")\nset_theme!(font = \"Arial\")\n\ngene = \"CACNA1G\"\nchr, start, stop = GeneticsMakie.findgene(gene, gencode)\nf = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs = GeneticsMakie.plotgenes!(ax, chr, start - 1e5, stop + 1e5, gencode)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 1e5, stop + 1e5)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"(Image: )","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"We can adjust the height of exons using the height keyword argument.","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"f = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs = GeneticsMakie.plotgenes!(ax, chr, start - 1e5, stop + 1e5, gencode; height = 0.1)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 1e5, stop + 1e5)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"(Image: )","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"We can change the color of genes and text using the genecolor and textcolor keyword arguments, respectively.","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"f = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs = GeneticsMakie.plotgenes!(ax, chr, start - 1e5, stop + 1e5, gencode; height = 0.1, genecolor = :mediumorchid3, textcolor = :forestgreen)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 1e5, stop + 1e5)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"(Image: )","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"Alternatively, we can visualize this locus by passing gene as a positional argument and  window as a keyword argument.","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"f = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs = GeneticsMakie.plotgenes!(ax, gene, gencode; window = 1e5, height = 0.1, genecolor = :brown3)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 1e5, stop + 1e5)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"(Image: )","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"There is no limit to the number of genes we can visualize. Below we visualize a larger  ± 2 Mb window around gene.","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"f = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs = GeneticsMakie.plotgenes!(ax, gene, gencode; window = 2e6, height = 0.1)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 2e6, stop + 2e6)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"(Image: )","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"note: Gene density\nAs some regions have higher gene density than the others, it would be wise (for publication purpose) to visualize a smaller genomic window for such gene-dense regions.","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"We can highlight a gene or sets of genes as below. This can be useful when highlighting genes  by certain properties such as those that are protein coding or those that are  loss-of-function intolerant or those that are  significant in some sort of gene-level association. ","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"f = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs = GeneticsMakie.plotgenes!(ax, chr, start, (gene, :brown3), gencode; window = 1e5, height = 0.1)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 1e5, start + 1e5)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"(Image: )","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"f = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs = GeneticsMakie.plotgenes!(ax, chr, start, ([gene, \"EPN3\"], [:brown3, :forestgreen]), gencode; window = 1e5, height = 0.1)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 1e5, start + 1e5)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"(Image: )","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"Finally, we can make additional modifications on top of the figure as needed using Makie.jl.","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"f = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs = GeneticsMakie.plotgenes!(ax, chr, start, gencode; window = 1e5, height = 0.1)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, start - 1e5, start + 1e5)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nvlines!(ax, start, color = (:gold, 0.5), linewidth = 0.5)\nvlines!(ax, stop, color = (:gold, 0.5), linewidth = 0.5)\nf","category":"page"},{"location":"examples/genes/","page":"Plotting genes","title":"Plotting genes","text":"(Image: )","category":"page"},{"location":"examples/locus/#Plotting-LocusZooom","page":"Plotting LocusZoom","title":"Plotting LocusZooom","text":"","category":"section"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"After Parsing GENCODE and Munging summary statistics, we can now put the pieces together to draw the backbone of a LocusZoom plot.  We can focus on JAZF1 locus as an example, which reaches strong genome-wide significance in GWAS for height. By default, GeneticsMakie.plotlocus! returns a straightforward scatter plot.","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"gene = \"JAZF1\"\nchr, start, stop = GeneticsMakie.findgene(gene, gencode)\nrange1 = start - 1e6\nrange2 = stop + 1e6\n\nn = length(df)\ntitles = [\"Height (Yengo et al. 2018)\", \"Weight (Yengo et al. 2018)\"]\nf = Figure(resolution = (306, 792))\naxs = [Axis(f[i, 1]) for i in 1:(n + 1)]\nfor i in 1:n\n    GeneticsMakie.plotlocus!(axs[i], chr, range1, range2, df[i])\n    rowsize!(f.layout, i, 30)\n    lines!(axs[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)\n    Label(f[i, 1, Top()], \"$(titles[i])\", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))\nend\nrs = GeneticsMakie.plotgenes!(axs[n + 1], chr, range1, range2, gencode; height = 0.1)\nrowsize!(f.layout, n + 1, rs)\nGeneticsMakie.labelgenome(f[n + 1, 1, Bottom()], chr, range1, range2)\nColorbar(f[1:n, 2], limits = (0, 1), ticks = 0:1:1, height = 20,\n    colormap = (:gray60, :red2), label = \"LD\", ticksize = 0, tickwidth = 0,\n    tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,\n    labelsize = 6, width = 5, spinewidth = 0.5)\nLabel(f[1:n, 0], text = \"-log[p]\", textsize = 6, rotation = pi / 2)\nrowgap!(f.layout, 5)\ncolgap!(f.layout, 5)\nfor i in 1:(n + 1)\n    vlines!(axs[i], start, color = (:gold, 0.5), linewidth = 0.5)\n    vlines!(axs[i], stop, color = (:gold, 0.5), linewidth = 0.5)\nend\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"(Image: )","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"To color variants by linkage disequilibrium (LD), we need a reference panel. If we already have one, we can use SnpArrays.jl to read in PLINK bed files. If not, we can download one as below. For the time being,  we only download and convert a single chromosome from the 1000 Genomes Project. ","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"using Pkg\nPkg.add(\"SnpArrays\")","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"using SnpArrays\n# Download 1000 Genomes data for a single chromosome\nbeagle = \"http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a\"\nurl = joinpath(beagle, \"b37.vcf/chr$(chr).1kg.phase3.v5a.vcf.gz\")\nvcf = basename(url)\nisdir(\"data/1kg\") || mkdir(\"data/1kg\")\nisfile(\"data/1kg/$(vcf)\") || download(url, \"data/1kg/$(vcf)\")\n# Convert vcf file to plink bed file (this step takes a while)\nisfile(\"data/1kg/$(replace(vcf, \".vcf.gz\" => \".bed\"))\") || vcf2plink(\"data/1kg/$(vcf)\", \"data/1kg/$(replace(vcf, \".vcf.gz\" => \"\"))\")\n# Download sample metadata\nurl = joinpath(beagle, \"sample_info/integrated_call_samples_v3.20130502.ALL.panel\")\nmeta = basename(url) \nisfile(\"data/1kg/$(meta)\") || download(url, \"data/1kg/$(meta)\")\n# Subset data to the genomic region of interest and European samples\nkgp = SnpData(\"data/1kg/$(replace(vcf, \".vcf.gz\" => \"\"))\")\nmeta = CSV.read(\"data/1kg/$(meta)\", DataFrame)\neur = meta.sample[meta.super_pop .== \"EUR\"]\ncolinds = findall((kgp.snp_info.position .>= range1) .& (kgp.snp_info.position .<= range2))\nrowinds = findall(in(eur), kgp.person_info.iid)\nfile = replace(vcf, \".vcf.gz\" => \".eur\")\nSnpArrays.filter(kgp, rowinds, colinds; des = \"data/1kg/$(file)\")\n# Apply minor allele frequency > 0.05 filter\nkgp = SnpData(\"data/1kg/$(file)\")\ncolinds = SnpArrays.filter(kgp.snparray; min_maf = 0.05)[2]\nfile = file * \".maf0.05\"\nSnpArrays.filter(kgp, trues(size(kgp)[1]), colinds; des = \"data/1kg/$(file)\")\nkgp = SnpData(\"data/1kg/$(file)\")","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"We can color variants by LD with the index/sentinel SNP by using the ld keyword argument.","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"f = Figure(resolution = (306, 792))\naxs = [Axis(f[i, 1]) for i in 1:(n + 1)]\nfor i in 1:n\n    GeneticsMakie.plotlocus!(axs[i], chr, range1, range2, df[i], ld = kgp)\n    rowsize!(f.layout, i, 30)\n    lines!(axs[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)\n    Label(f[i, 1, Top()], \"$(titles[i])\", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))\nend\nrs = GeneticsMakie.plotgenes!(axs[n + 1], chr, range1, range2, gencode; height = 0.1)\nrowsize!(f.layout, n + 1, rs)\nGeneticsMakie.labelgenome(f[n + 1, 1, Bottom()], chr, range1, range2)\nColorbar(f[1:n, 2], limits = (0, 1), ticks = 0:1:1, height = 20,\n    colormap = (:gray60, :red2), label = \"LD\", ticksize = 0, tickwidth = 0,\n    tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,\n    labelsize = 6, width = 5, spinewidth = 0.5)\nLabel(f[1:n, 0], text = \"-log[p]\", textsize = 6, rotation = pi / 2)\nrowgap!(f.layout, 5)\ncolgap!(f.layout, 5)\nfor i in 1:(n + 1)\n    vlines!(axs[i], start, color = (:gold, 0.5), linewidth = 0.5)\n    vlines!(axs[i], stop, color = (:gold, 0.5), linewidth = 0.5)\nend\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"(Image: )","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"We can also color variants by LD with the same SNP by using the ld keyword argument.","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"f = Figure(resolution = (306, 792))\naxs = [Axis(f[i, 1]) for i in 1:(n + 1)]\nfor i in 1:n\n    GeneticsMakie.plotlocus!(axs[i], chr, range1, range2, df[i], ld = (kgp, \"rs508347\"))\n    rowsize!(f.layout, i, 30)\n    lines!(axs[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)\n    Label(f[i, 1, Top()], \"$(titles[i])\", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))\nend\nrs = GeneticsMakie.plotgenes!(axs[n + 1], chr, range1, range2, gencode; height = 0.1)\nrowsize!(f.layout, n + 1, rs)\nGeneticsMakie.labelgenome(f[n + 1, 1, Bottom()], chr, range1, range2)\nColorbar(f[1:n, 2], limits = (0, 1), ticks = 0:1:1, height = 20,\n    colormap = (:gray60, :red2), label = \"LD\", ticksize = 0, tickwidth = 0,\n    tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,\n    labelsize = 6, width = 5, spinewidth = 0.5)\nLabel(f[1:n, 0], text = \"-log[p]\", textsize = 6, rotation = pi / 2)\nrowgap!(f.layout, 5)\ncolgap!(f.layout, 5)\nfor i in 1:(n + 1)\n    vlines!(axs[i], start, color = (:gold, 0.5), linewidth = 0.5)\n    vlines!(axs[i], stop, color = (:gold, 0.5), linewidth = 0.5)\nend\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"(Image: )","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"By using Makie.jl's layout tools,  it becomes easy to draw additional tracks. For example, in a separate track,  the variants could be colored or could have varying sizes depending on their minor allele frequency.  In another example, the variants could be colored based on their inclusion in a  credible set post-fine-mapping.","category":"page"},{"location":"examples/locus/","page":"Plotting LocusZoom","title":"Plotting LocusZoom","text":"note: Plotting the intersection of SNPs, not the union\nGeneticsMakie.plotlocus! plots only the variants that are present in the reference panel,  when the ld keyword argument is specified. Although SNPs that are missing in the reference panel could be plotted differently (e.g. with varying transparency and shape), GeneticsMakie.jl is designed to visualize 100s of phenotypes simultaneously in which case such discrepancy is hard to tell and  is confusing. Hence, for more direct comparison of loci across phenotypes,  only the variants that are found in the reference panel are shown.","category":"page"},{"location":"examples/gtf/#Parsing-GENCODE","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"","category":"section"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"Install the relevant packages in the usual way.","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"using Pkg\nPkg.add([\"GeneticsMakie\", \"CSV\", \"DataFrames\"])","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"To plot genes and isoforms, we need a transcriptome annotation. We can use  the latest GENCODE annotation for  the human genome (GRCh37), where we download the comprehensive  gene annotation file in GTF format.","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"using GeneticsMakie, CSV, DataFrames\nisdir(\"data\") || mkdir(\"data\")\nurl = \"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/gencode.v39lift37.annotation.gtf.gz\"\nfile = basename(url)\nisdir(\"data/gencode\") || mkdir(\"data/gencode\")\nisfile(\"data/gencode/$(file)\") || download(url, \"data/gencode/$(file)\")\nh = [\"seqnames\", \"source\", \"feature\", \"start\", \"end\", \"score\", \"strand\", \"phase\", \"info\"]\ngencode = CSV.read(\"data/gencode/$(file)\", DataFrame; delim = \"\\t\", comment = \"#\", header = h)","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"warning: Human genome build\nThe latest human genome assembly is GRCh38, but we use an annotation with coordinates  from the older version (GRCh37), because a lot of the GWAS results are shared in  GRCh37 genomic coordinates. Make sure to use the matching human genome build when visualizing your results. ","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"The ninth column of a GTF file  contains rich information about features, so we can parse this column.","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"GeneticsMakie.parsegtf!(gencode)","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"info: Chromosome names\nChromosome names are munged to not contain “chr” prefix, and their type is String, since there could be non-numerical chromosome names, such as sex chromosomes and mitochondrial genome.","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"To reduce memory intake, we can also subset gencode to most commonly used columns in downstream analyses.","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"select!(gencode, :seqnames, :feature, :start, :end, :strand, :gene_name, :gene_type, :transcript_id)","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"Other transcriptome annotations, such as one from RefSeq, can be used for plotting functions  as long as they contain the above columns with the right column names.","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"Once gencode is ready, we can look up where the gene is on the human genome.","category":"page"},{"location":"examples/gtf/","page":"Parsing GENCODE","title":"Parsing GENCODE","text":"GeneticsMakie.findgene(\"RBFOX1\", gencode)","category":"page"},{"location":"examples/isoforms/#Plotting-isoforms","page":"Plotting isoforms","title":"Plotting isoforms","text":"","category":"section"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"We can focus on NRXN1 gene as our initial example. GeneticsMakie.plotisoforms! returns  genomic coordinates for the gene of interest so that an appropriate label can be passed onto  GeneticsMakie.labelgenome. NRXN1 gene has many isoforms as we see below, and even more  isoforms are likely to be discovered in the future. For this reason, plotting isoforms of multiple genes is not available. ","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"gene = \"NRXN1\"\nf = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs, chr, range1, range2 = GeneticsMakie.plotisoforms!(ax, gene, gencode; height = 0.1)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, range1, range2)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"(Image: )","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"To save some space, we can plot the isoform labels on the left by using the text keyword argument. ","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"gene = \"NRXN1\"\nf = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs, chr, range1, range2 = GeneticsMakie.plotisoforms!(ax, gene, gencode; height = 0.1, text = :l)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, range1, range2)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"(Image: )","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"We can change the color of isoforms and text using the isoformcolor and textcolor keyword arguments, respectively. We look at a different gene, GRIN2A.","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"gene = \"GRIN2A\"\nf = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs, chr, range1, range2 = GeneticsMakie.plotisoforms!(ax, gene, gencode; isoformcolor = :forestgreen, height = 0.1)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, range1, range2)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"(Image: )","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"We can change the order of isoforms by using the orderby keyword argument.","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"f = Figure(resolution = (306, 792))\nax = Axis(f[1, 1])\nrs, chr, range1, range2 = GeneticsMakie.plotisoforms!(ax, gene, gencode; orderby = [\"ENST00000675189\", \"ENST00000463531\"], isoformcolor = :forestgreen, height = 0.1)\nGeneticsMakie.labelgenome(f[1, 1, Bottom()], chr, range1, range2)\nrowsize!(f.layout, 1, rs)\nresize_to_layout!(f)\nf","category":"page"},{"location":"examples/isoforms/","page":"Plotting isoforms","title":"Plotting isoforms","text":"(Image: )","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GeneticsMakie","category":"page"},{"location":"#GeneticsMakie","page":"Home","title":"GeneticsMakie","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The goal of GeneticsMakie.jl is to  permit seamless data visualization and exploratory data analysis of the human genome within the larger Julia data science and OpenMendel ecosystems. The package provides convenient wrapper functions for wrangling genetic association results and  plotting them using Makie.jl. Every component of a figure  can be easily customized and extended, and the package generates high-quality, publication-ready figures. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: \"mhc\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"tip: Getting started\nPlease peruse the documentations of  Makie.jl,  CSV.jl,  DataFrames.jl,  and SnpArrays.jl.  Familiarity with these packages will allow visualization of most types of genetic and genomic data.  Makie.jl's default layout tools are particularly useful for  plotting different genetic and genomic data modalities as separate layers.","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: An usage case\nIf you have run a genome-wide association study (GWAS) at the variant-level,  and you would like to eyeball genome-wide significant loci across hundreds of phenotypes, then you are in the right place.","category":"page"},{"location":"examples/summary/#Munging-summary-statistics","page":"Munging summary statistics","title":"Munging summary statistics","text":"","category":"section"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"GWAS summary statistics come in a variety of shapes and flavors, so harmonizing them is crucial in making our lives easier when trying to visualize their results. We can  take a peak at GWAS results for height and weight, the two classic anthropometric traits. ","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"using Pkg\nPkg.add([\"GeneticsMakie\", \"CSV\", \"DataFrames\"])","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"using GeneticsMakie, CSV, DataFrames\ngwas = Dict(\n    \"height\" => \"https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz\",\n    \"weight\" => \"https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz\"\n)\nisdir(\"data/gwas\") || mkdir(\"data/gwas\")\ndf = DataFrame[]\nfor key in keys(gwas)\n    url = gwas[key]\n    isfile(\"data/gwas/$(basename(url))\") || download(url, \"data/gwas/$(basename(url))\")\n    push!(df, CSV.read(\"data/gwas/$(basename(url))\", DataFrame; comment = \"##\", missingstring = [\"NA\"]))\nend","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"To harmonize summary statistics, we run a single command GeneticsMakie.mungesumstats!.  For all downstream plotting functions, we note that GWAS summary statistics should be pre-harmonized.","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"GeneticsMakie.mungesumstats!(df)","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"To find GWAS loci for each phenotype that are separated by at least 1 Mb,  we can execute the following command.","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"GeneticsMakie.findgwasloci(df[1])\nGeneticsMakie.findgwasloci(df[2])","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"To find every GWAS loci across multiple phenotypes that are separated by at least 1 Mb, we can instead run the following command.","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"GeneticsMakie.findgwasloci(df)","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"Such an exhaustive list of loci can then be iterated through and visualized by  Plotting LocusZooom.","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"We can also find the closest gene to each index SNP in GWAS loci.","category":"page"},{"location":"examples/summary/","page":"Munging summary statistics","title":"Munging summary statistics","text":"loci = GeneticsMakie.findgwasloci(df[1])\nGeneticsMakie.findclosestgene(loci, gencode)\nGeneticsMakie.findclosestgene(loci, gencode; start = true) # closest gene from gene start site\nGeneticsMakie.findclosestgene(loci, gencode; proteincoding = true) # closest \"protein-coding\" gene","category":"page"}]
}
