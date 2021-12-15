"""
    coordinategenes(chromosome::AbstractString, range1::Real, range2::Real, gencode::DataFrame)

Subset `gencode` to a given `chromosome` and genomic range between `range1` and
`range2`, and determine coordinates of exons for each gene residing in the genomic region.
"""
function coordinategenes(chromosome::AbstractString, 
    range1::Real, 
    range2::Real, 
    gencode::DataFrame,
    height::Real)

    df = filter(x -> (x.seqnames == chromosome) && (x.end >= range1) && (x.start <= range2), gencode)
    dfg = df[df.feature .== "gene", :]
    dfe = df[df.feature .== "exon", :]
    genes = unique(dfg.gene_name)
    strand = [dfg.strand[findfirst(isequal(gene), dfg.gene_name)] for gene in genes]
    ps = Vector{Vector{Polygon}}(undef, length(genes))
    bs = Matrix{Float64}(undef, length(genes), 2)
    rows = ones(Int64, length(genes))
    for j in eachindex(genes)
        ind = findfirst(isequal(genes[j]), dfg.gene_name)
        start1 = dfg.start[ind]
        stop1 = dfg[ind, :end]
        center1 = (start1 + stop1) / 2
        label1 = center1 + (length(genes[j]) + 1) / 2 * 44000
        bs[j, 1] = start1
        bs[j, 2] = stop1
        for k in (j + 1):length(genes)
            ind = findfirst(isequal(genes[k]), dfg.gene_name)
            start2 = dfg.start[ind]
            stop2 = dfg[ind, :end]
            center2 = (start2 + stop2) / 2
            label2 = center2 - (length(genes[k]) + 1)/ 2 * 44000
            if ((stop1 > start2) || (label1 > start2) || (label1 > label2)) && (rows[j] == rows[k])
                rows[k] = rows[j] + 1
                for l in 1:(j - 1)
                    start3 = bs[l, 1]
                    stop3 = bs[l, 2]
                    center3 = (start3 + stop3) / 2
                    label3 = center3 + (length(genes[l]) + 1) / 2 * 44000
                    if ((stop3 > start2) || (label3 > start2) || (label3 > label2)) && (rows[l] == rows[k])
                        rows[k] = rows[l] + 1
                    end
                end
            end
        end
    end
    for j in eachindex(genes)
        ranges = Matrix{Float32}(dfe[findall(isequal(genes[j]), dfe.gene_name), [:start, :end]])
        n = size(ranges, 1)
        p = Vector{Polygon}(undef, n)
        for i = 1:n
            p[i] = Polygon(
                [Point2f(ranges[i, 1], 1 - height - (rows[j] - 1) * (0.25 + height)),
                Point2f(ranges[i, 1], 1 - (rows[j] - 1) * (0.25 + height)),
                Point2f(ranges[i, 2], 1 - (rows[j] - 1) * (0.25 + height)),
                Point2f(ranges[i, 2], 1 - height - (rows[j] - 1) * (0.25 + height))]
            )
        end
        ps[j] = p
    end
    return genes, strand, ps, bs, rows
end

"""
    plotgenes!(ax::Axis, chromosome::AbstractString, range1::Real, range2::Real, gencode::DataFrame; height, highlight)

Plot collapsed gene bodies for genes within a given `chromosome` and genomic range 
between `range1` and `range2`. Optionally, height of exons can be adjusted using
`height` and genes can be highlighted via `highlight`.
"""
function plotgenes!(ax::Axis,
    chromosome::AbstractString, 
    range1::Real,
    range2::Real,
    gencode::DataFrame;
    height::Real = 0.25,
    highlight::Union{Nothing, AbstractVector{<:AbstractString}} = nothing)

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
    if isnothing(highlight)
        for j in 1:size(ps, 1)
            poly!(ax, ps[j], color = :royalblue, strokewidth = 0)
            lines!(ax, [bs[j, 1], bs[j, 2]], 
                [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)], 
                color = :royalblue, linewidth = 0.5)
            g = (strand[j] == "+" ? genes[j] * "→" : "←" * genes[j])
            text!(ax, "$g", 
                position = ((bs[j, 1] + bs[j, 2]) / 2, 1 - (rows[j] - 1) * (0.25 + height)), 
                align = (:center, :bottom), textsize = 6)
        end
    else
        for j in 1:size(ps, 1)
            if genes[j] in highlight
                poly!(ax, ps[j], color = "#3484BD", strokewidth = 0)
                lines!(ax, [bs[j, 1], bs[j, 2]], 
                    [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)], 
                    color = "#3484BD", linewidth = 0.5)
                g = (strand[j] == "+" ? genes[j] * "→" : "←" * genes[j])
                text!(ax, "$g", 
                    position = ((bs[j, 1] + bs[j, 2]) / 2, 1 - (rows[j] - 1) * (0.25 + height)), 
                    align = (:center, :bottom), textsize = 6, color = "#3484BD")
            else
                poly!(ax, ps[j], color = "#BEBEBE", strokewidth = 0)
                lines!(ax, [bs[j, 1], bs[j, 2]], 
                    [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)], 
                    color = "#BEBEBE", linewidth = 0.5)
                g = (strand[j] == "+" ? genes[j] * "→" : "←" * genes[j])
                text!(ax, "$g", 
                    position = ((bs[j, 1] + bs[j, 2]) / 2, 1 - (rows[j] - 1) * (0.25 + height)),
                    align = (:center, :bottom), textsize = 6, color = "#BEBEBE")
            end
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


"""
    plotgenes!(ax::Axis, chromosome::AbstractString, bp::Real, gencode::DataFrame; window, height, highlight)

Plot collapsed gene bodies for genes within a given `chromosome` and a certain
`window` around a genomic coordinate `bp`. The default window is 1 Mb. Optionally, 
height of exons can be adjusted using `height` and genes can be highlighted via `highlight`.
"""
plotgenes!(ax::Axis, chromosome::AbstractString, bp::Real, gencode::DataFrame; window::Real = 1e6, kwargs...) =
    plotgenes!(ax, chromosome, bp - window, bp + window, gencode; kwargs...)

"""
    plotgenes!(ax::Axis, gene::AbstractString, gencode::DataFrame; window, height, highlight)

Plot collapsed gene bodies for genes within a certain window around `gene`. 
The default window is 1 Mb. Optionally, height of exons can be adjusted using 
`height` and genes can be highlighted via `highlight`.
"""
function plotgenes!(ax::Axis, gene::AbstractString, gencode::DataFrame; window::Real = 1e6, kwargs...)
    ind = findfirst(isequal(gene), gencode.gene_name)
    plotgenes!(ax, gencode.seqnames[ind], gencode.start[ind] - window, gencode[ind, :end] + window, gencode; kwargs...)
end