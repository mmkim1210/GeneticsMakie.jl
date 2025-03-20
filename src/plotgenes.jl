function coordinategenes(
    chromosome::AbstractString,
    range1::Real,
    range2::Real,
    gencode::DataFrame,
    height::Real
    )
    df = filter(x -> (x.seqnames == chromosome) && (x.end >= range1) && (x.start <= range2), gencode)
    dfg = view(df, df.feature .== "gene", :)
    dfe = view(df, df.feature .== "exon", :)
    genes = unique(dfg.gene_name)
    strand = [dfg.strand[findfirst(isequal(gene), dfg.gene_name)] for gene in genes]
    ps = Vector{Vector{Rect}}(undef, length(genes))
    bs = Matrix{Float64}(undef, length(genes), 2)
    rows = ones(Int, length(genes))
    prop =  23000 * (range2 - range1) / 2.2e6
    for j in eachindex(genes)
        ind = findfirst(isequal(genes[j]), dfg.gene_name)
        start1 = dfg.start[ind]
        stop1 = dfg[ind, :end]
        center1 = (start1 + stop1) / 2
        label1 = center1 + (length(genes[j]) + 1) * prop
        bs[j, 1] = start1
        bs[j, 2] = stop1
        for k in (j + 1):length(genes)
            ind = findfirst(isequal(genes[k]), dfg.gene_name)
            start2 = dfg.start[ind]
            stop2 = dfg[ind, :end]
            center2 = (start2 + stop2) / 2
            label2 = center2 - (length(genes[k]) + 1) * prop
            if ((stop1 > start2) || (label1 > label2)) && (rows[j] == rows[k])
                rows[k] = rows[j] + 1
                while true
                    l = findprev(isequal(rows[k]), rows, k - 1)
                    isnothing(l) && break
                    start3 = bs[l, 1]
                    stop3 = bs[l, 2]
                    center3 = (start3 + stop3) / 2
                    label3 = center3 + (length(genes[l]) + 1) * prop
                    ((stop3 < start2) && (label3 < label2)) && break
                    rows[k] = rows[l] + 1
                end
            end
        end
    end
    for j in eachindex(genes)
        ranges = view(dfe, findall(isequal(genes[j]), dfe.gene_name), [:start, :end])
        n = size(ranges, 1)
        p = Vector{Rect}(undef, n)
        for i = 1:n
            p[i] = Rect([
                Point2f(ranges[i, 1], 1 - height - (rows[j] - 1) * (0.25 + height)),
                Point2f(ranges[i, 1], 1 - (rows[j] - 1) * (0.25 + height)),
                Point2f(ranges[i, 2], 1 - (rows[j] - 1) * (0.25 + height)),
                Point2f(ranges[i, 2], 1 - height - (rows[j] - 1) * (0.25 + height))
            ])
        end
        ps[j] = p
    end
    return genes, strand, ps, bs, rows
end

"""
    plotgenes!(ax::Axis, chromosome::AbstractString, range1::Real, range2::Real, gencode::DataFrame)
    plotgenes!(ax::Axis, chromosome::AbstractString, bp::Real, gencode::DataFrame)
    plotgenes!(ax::Axis, gene::AbstractString, gencode::DataFrame)

Plot collapsed gene bodies for genes within a given `chromosome` and genomic range
between `range1` and `range2`.

Alternatively, plot within a given `chromosome` and a certain `window` around a
genomic coordinate `bp` or plot within a certain `window` around `gene`.

# Keyword arguments
```
height      height of exons; default 0.25
genecolor   color of genes; default :royalblue
textcolor   color of gene labels; default :black
window      window around bp or gene; default 1e6
```
"""
function plotgenes!(
    ax::Axis,
    chromosome::AbstractString,
    range1::Real,
    range2::Real,
    gencode::DataFrame;
    height::Real = 0.25,
    genecolor = :royalblue,
    textcolor = :black
    )
    genes, strand, ps, bs, rows = coordinategenes(chromosome, range1, range2, gencode, height)
    if length(rows) == 0
        ax.spinewidth = 0.75
        hidexdecorations!(ax)
        hideydecorations!(ax)
        xlims!(ax, range1, range2)
        ylims!(ax, 0.875 - height, 1.375)
        rs = 18 * (0.25 + height) / 0.5
        return rs
    end
    for j in 1:size(ps, 1)
        poly!(ax, ps[j], color = genecolor, strokewidth = 0)
        lines!(ax, [bs[j, 1], bs[j, 2]],
            [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)],
            color = genecolor, linewidth = 0.5)
        g = (strand[j] == "+" ? genes[j] * "→" : "←" * genes[j])
        text!(ax, (bs[j, 1] + bs[j, 2]) / 2, 1 - (rows[j] - 1) * (0.25 + height),
            text = "$g", align = (:center, :bottom), fontsize = 6, color = textcolor)
    end
    ax.spinewidth = 0.75
    hidexdecorations!(ax)
    hideydecorations!(ax)
    xlims!(ax, range1, range2)
    ylims!(ax, 0.875 - height  - (maximum(rows) - 1) * (0.25 + height), 1.375)
    rs = 18 * maximum(rows) * (0.25 + height) / 0.5
    return rs
end

plotgenes!(ax::Axis, chromosome::AbstractString, bp::Real, gencode::DataFrame; window::Real = 1e6, kwargs...) =
    plotgenes!(ax, chromosome, bp - window, bp + window, gencode; kwargs...)

function plotgenes!(ax::Axis, gene::AbstractString, gencode::DataFrame; window::Real = 1e6, kwargs...)
    chr, start, stop = findgene(gene::AbstractString, gencode::DataFrame)
    plotgenes!(ax, chr, start - window, stop + window, gencode; kwargs...)
end

"""
    plotgenes!(ax::Axis, chromosome::AbstractString, range1::Real, range2::Real, highlight::Tuple{AbstractVector, AbstractVector}, gencode::DataFrame; height::Real)
    plotgenes!(ax::Axis, chromosome::AbstractString, bp::Real, highlight::Tuple{AbstractVector, AbstractVector}, gencode::DataFrame; window::Real, height::Real)
    plotgenes!(ax::Axis, gene::AbstractString, highlight::Tuple{AbstractVector, AbstractVector}, gencode::DataFrame; window::Real, height::Real)

Plot gene bodies with a vector of genes highlighted by a vector of colors via `highlight`.
"""
function plotgenes!(
    ax::Axis,
    chromosome::AbstractString,
    range1::Real,
    range2::Real,
    highlight::Tuple{AbstractVector, AbstractVector},
    gencode::DataFrame;
    height::Real = 0.25
)

    genes, strand, ps, bs, rows = coordinategenes(chromosome, range1, range2, gencode, height)
    if length(rows) == 0
        ax.spinewidth = 0.75
        hidexdecorations!(ax)
        hideydecorations!(ax)
        xlims!(ax, range1, range2)
        ylims!(ax, 0.875 - height, 1.375)
        rs = 18 * (0.25 + height) / 0.5
        return rs
    end
    for j in 1:size(ps, 1)
        ind = findfirst(isequal(genes[j]), highlight[1])
        if !isnothing(ind)
            poly!(ax, ps[j], color = highlight[2][ind], strokewidth = 0)
            lines!(ax, [bs[j, 1], bs[j, 2]],
                [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)],
                color = highlight[2][ind], linewidth = 0.5)
            g = (strand[j] == "+" ? genes[j] * "→" : "←" * genes[j])
            text!(ax, (bs[j, 1] + bs[j, 2]) / 2, 1 - (rows[j] - 1) * (0.25 + height),
                text = "$g", align = (:center, :bottom), fontsize = 6, color = highlight[2][ind])
        else
            poly!(ax, ps[j], color = "gray60", strokewidth = 0)
            lines!(ax, [bs[j, 1], bs[j, 2]],
                [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)],
                color = "gray60", linewidth = 0.5)
            g = (strand[j] == "+" ? genes[j] * "→" : "←" * genes[j])
            text!(ax, (bs[j, 1] + bs[j, 2]) / 2, 1 - (rows[j] - 1) * (0.25 + height),
                text = "$g", align = (:center, :bottom), fontsize = 6, color = "gray60")
        end
    end
    ax.spinewidth = 0.75
    hidexdecorations!(ax)
    hideydecorations!(ax)
    xlims!(ax, range1, range2)
    ylims!(ax, 0.875 - height  - (maximum(rows) - 1) * (0.25 + height), 1.375)
    rs = 18 * maximum(rows) * (0.25 + height) / 0.5
    return rs
end

plotgenes!(ax::Axis, chromosome::AbstractString, bp::Real, highlight::Tuple{AbstractVector, AbstractVector}, gencode::DataFrame; window::Real = 1e6, kwargs...) =
    plotgenes!(ax, chromosome, bp - window, bp + window, highlight, gencode; kwargs...)

function plotgenes!(ax::Axis, gene::AbstractString, highlight::Tuple{AbstractVector, AbstractVector}, gencode::DataFrame; window::Real = 1e6, kwargs...)
    chr, start, stop = findgene(gene::AbstractString, gencode::DataFrame)
    plotgenes!(ax, chr, start - window, stop + window, highlight, gencode; kwargs...)
end

plotgenes!(ax::Axis, chromosome::AbstractString, range1::Real, range2::Real, highlight::Tuple{AbstractString, Any}, gencode::DataFrame; kwargs...) =
    plotgenes!(ax, chromosome, range1, range2, ([highlight[1]], [highlight[2]]), gencode; kwargs...)

plotgenes!(ax::Axis, chromosome::AbstractString, bp::Real, highlight::Tuple{AbstractString, Any}, gencode::DataFrame; kwargs...) =
    plotgenes!(ax, chromosome, bp, ([highlight[1]], [highlight[2]]), gencode; kwargs...)

plotgenes!(ax::Axis, gene::AbstractString, highlight::Tuple{AbstractString, Any}, gencode::DataFrame; kwargs...) =
    plotgenes!(ax, gene, ([highlight[1]], [highlight[2]]), gencode; kwargs...)

"""
    labelgenome(g::GridPosition, chromosome::AbstractString, range1::Real, range2::Real)

Label `g` with a given `chromosome` and genomic range between `range1` and `range2`.
"""
function labelgenome(g::GridPosition, chromosome::AbstractString, range1::Real, range2::Real)
    Label(g, "~$(round(range1 / 1e6; digits = 1)) Mb", fontsize = 6, halign = :left)
    Label(g, "Chr $(chromosome)", fontsize = 6, halign = :center)
    Label(g, "~$(round(range2 / 1e6; digits = 1)) Mb", fontsize = 6, halign = :right)
end
