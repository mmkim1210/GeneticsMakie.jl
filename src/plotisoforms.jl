"""
    coordinateisforms(gene::AbstractString, gencode::DataFrame, orderby::Union{Nothing, AbstractVector{<:AbstractString}}, height::Real, text::Bool, textpos::Symbol)

Subset `gencode` to a `gene` and determine coordinates of exons for each isoform.
"""
function coordinateisforms(gene::AbstractString, 
    gencode::DataFrame,
    orderby::Union{Nothing, AbstractVector{<:AbstractString}},
    height::Real,
    text::Bool,
    textpos::Symbol)

    df = filter(x -> x.gene_name == gene, gencode)
    dfi = view(df, df.feature .== "transcript", :)
    dfe = view(df, df.feature .== "exon", :)
    chromosome = df.seqnames[1]
    isoforms = unique(dfi.transcript_id)
    n = length(isoforms)
    ps = Vector{Vector{Polygon}}(undef, n) 
    bs = Matrix{Float64}(undef, n, 2)
    rows = collect(1:n)
    if !isnothing(orderby)
        dforder = DataFrame(transcript_id = [orderby; filter(!in(orderby), isoforms)], rank = 1:n)
        dfi = leftjoin!(dfi, dforder; on = :transcript_id)
        sort!(dfi, [order(:rank)])
        isoforms = unique(dfi.transcript_id)
    end
    for j in eachindex(isoforms)
        ranges = view(dfe, findall(isequal(isoforms[j]), dfe.transcript_id), [:start, :end])
        m = size(ranges, 1)
        p = Vector{Polygon}(undef, m)
        if text && textpos == :top
            for i = 1:m
                p[i] = Polygon(
                    [Point2f(ranges[i, 1], 1 - height - (rows[j] - 1) * (0.25 + height)),
                    Point2f(ranges[i, 1], 1 - (rows[j] - 1) * (0.25 + height)),
                    Point2f(ranges[i, 2], 1 - (rows[j] - 1) * (0.25 + height)),
                    Point2f(ranges[i, 2], 1 - height - (rows[j] - 1) * (0.25 + height))]
                )
            end
        else
            for i = 1:m
                p[i] = Polygon(
                    [Point2f(ranges[i, 1], 1 - height - (rows[j] - 1) * (0.025 + height)),
                    Point2f(ranges[i, 1], 1 - (rows[j] - 1) * (0.025 + height)),
                    Point2f(ranges[i, 2], 1 - (rows[j] - 1) * (0.025 + height)),
                    Point2f(ranges[i, 2], 1 - height - (rows[j] - 1) * (0.025 + height))]
                )
            end
        end
        ps[j] = p
        ind = findfirst(isequal(isoforms[j]), dfi.transcript_id)
        start = dfi.start[ind]
        stop = dfi[ind, :end]
        bs[j, 1] = start
        bs[j, 2] = stop
    end
    return isoforms, ps, bs, rows, chromosome
end

"""
    plotisoforms!(ax::Axis, gene::AbstractString, gencode::DataFrame; orderby, height, text, textpos)

Plot each isoform of a given `gene` in a separate row. Optionally, order of
isoforms can be changed by `orderby`, height of exons can be adjusted using
`height`, and color of exons or gene names adjusted using `exoncolor` and `textcolor`,
respectively. The presence or absence of text label of isoforms can be specified
using `text` with its position `textpos`.
"""
function plotisoforms!(ax::Axis,
    gene::AbstractString,
    gencode::DataFrame;
    orderby::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    height::Real = 0.25,
    exoncolor::Union{Symbol, AbstractString} = :royalblue,
    textcolor::Union{Symbol, AbstractString} = :black,
    text::Bool = true,
    textpos::Symbol = :top)

    isoforms, ps, bs, rows, chromosome = coordinateisforms(gene, gencode, orderby, height, text, textpos)
    if text && textpos == :top
        for j in 1:size(ps, 1)
            poly!(ax, ps[j], color = exoncolor, strokewidth = 0)
            lines!(ax, [bs[j, 1], bs[j, 2]], 
                [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)],
                color = exoncolor, linewidth = 0.5)
            text!(ax, "$(isoforms[j])", 
                position = ((bs[j, 1] + bs[j, 2]) / 2, 1 - (rows[j] - 1) * (0.25 + height)), 
                align = (:center, :bottom), textsize = 6, color = textcolor)
        end
    else
        for j in 1:size(ps, 1)
            poly!(ax, ps[j], color = exoncolor, strokewidth = 0)
            lines!(ax, [bs[j, 1], bs[j, 2]], 
                [1 - height / 2 - (rows[j] - 1) * (0.025 + height), 1 - height / 2 - (rows[j] - 1) * (0.025 + height)],
                color = exoncolor, linewidth = 0.5)
        end
    end
    range1 = minimum(bs[:, 1])
    range2 = maximum(bs[:, 2])
    diff = range2 - range1
    prop = 23000 * diff / 2.2e6
    if text && textpos == :top
        storage = bs[:, 1] + bs[:, 2]
        ind = argmin(storage)
        range1 = storage[ind] / 2 - prop * 12
        ind = argmax(storage)
        range2 = storage[ind] / 2 + prop * 12
    else 
        range1 = range1 - diff / 50
        range2 = range2 + diff / 50
    end
    ax.spinewidth = 0.75
    hidexdecorations!(ax)
    if text && textpos != :top
        ax.yticks = ([1 - height / 2 - (rows[j] - 1) * (0.025 + height) for j in 1:length(isoforms)], isoforms)
        hideydecorations!(ax, ticklabels = false)
        ax.yticklabelsize = 4
        ax.yticklabelpad = 3
        textpos == :right ? ax.yaxisposition = :right : ax.yaxisposition = :left
        rs = 18 * maximum(rows) * (0.025 + height) / 0.5
    elseif text && textpos == :top
        hideydecorations!(ax)
        rs = 18 * maximum(rows) * (0.25 + height) / 0.5
    else
        hideydecorations!(ax)
        rs = 18 * maximum(rows) * (0.025 + height) / 0.5
    end
    xlims!(ax, range1, range2)
    if text && textpos == :top
        ylims!(ax, 0.875 - height  - (maximum(rows) - 1) * (0.25 + height), 1.375)
    else
        ylims!(ax, 0.975 - height  - (maximum(rows) - 1) * (0.025 + height), 1.05)
    end
    return rs, chromosome, range1, range2
end