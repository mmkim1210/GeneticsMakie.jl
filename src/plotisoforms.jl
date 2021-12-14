"""
    coordinateisforms(gene::AbstractString, gencode::DataFrame)

Subset `gencode` to a `gene` and determine coordinates of exons for each isoform.
"""
function coordinateisforms(gene::AbstractString, 
    gencode::DataFrame, 
    orderby::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    height::Real = 0.25,
    label::Bool = true,
    labelpos::Symbol = :top)

    df = filter(x -> x.gene_name == gene, gencode)
    dfi = df[df.feature .== "transcript", :]
    dfe = df[df.feature .== "exon", :]
    chromosome = df.seqnames[1]
    isoforms = unique(dfi.transcript_id)
    n = length(isoforms)
    ps = Vector{Vector{Polygon}}(undef, n) 
    bs = Matrix{Float64}(undef, n, 2)
    rows = collect(1:n)
    if !isnothing(orderby)
        dforder = DataFrame(transcript_id = [orderby; filter(!in(orderby), isoforms)], rank = 1:n)
        dfi = @chain dfi begin
            leftjoin(_, dforder; on = :transcript_id)
        end
        sort!(dfi, [order(:rank)])
        isoforms = unique(dfi.transcript_id)
    end
    for j in eachindex(isoforms)
        ranges = Matrix{Float32}(dfe[findall(isequal(isoforms[j]), dfe.transcript_id), [:start, :end]])
        m = size(ranges, 1)
        p = Vector{Polygon}(undef, m)
        if label && labelpos == :top
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
    plotisoforms(gene::AbstractString, gencode::DataFrame)

Plot each isoform of a given `gene` in a separate row.
"""
function plotisoforms(gene::AbstractString,
    gencode::DataFrame;
    filename::AbstractString = "isoforms",
    orderby::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    height::Real = 0.25,
    label::Bool = true,
    labelpos::Symbol = :top)

    isoforms, ps, bs, rows, chromosome = coordinateisforms(gene, gencode, orderby, height, label, labelpos)
    CairoMakie.activate!(type = "pdf")
    set_theme!(font = "Arial")
    f = Figure(resolution = (306, 792))
    ga = f[1, 1] = GridLayout()
    gb = f[2, 1] = GridLayout()
    ax = Axis(ga[1, 1])
    if label && labelpos == :top
        for j in 1:size(ps, 1)
            poly!(ax, ps[j], color = :royalblue, strokewidth = 0)
            lines!(ax, [bs[j, 1], bs[j, 2]], 
                [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)], 
                color = :royalblue, linewidth = 0.5)
            text!(ax, "$(isoforms[j])", 
                position = ((bs[j, 1] + bs[j, 2]) / 2, 1 - (rows[j] - 1) * (0.25 + height)), 
                align = (:center, :bottom), textsize = 6)
        end
    else
        for j in 1:size(ps, 1)
            poly!(ax, ps[j], color = :royalblue, strokewidth = 0)
            lines!(ax, [bs[j, 1], bs[j, 2]], 
                [1 - height / 2 - (rows[j] - 1) * (0.025 + height), 1 - height / 2 - (rows[j] - 1) * (0.025 + height)], 
                color = :royalblue, linewidth = 0.5)
        end
    end
    range1 = minimum(bs[:, 1])
    range2 = maximum(bs[:, 2])
    diff = range2 - range1
    range1 = range1 - diff / 20 
    range2 = range2 + diff / 20 
    ax.spinewidth = 0.75
    hidexdecorations!(ax)
    if label && labelpos != :top
        ax.yticks = ([1 - height / 2 - (rows[j] - 1) * (0.025 + height) for j in 1:length(isoforms)], isoforms)
        hideydecorations!(ax, ticklabels = false)
        ax.yticklabelsize = 4
        ax.yticklabelpad = 3
    else
        hideydecorations!(ax)
    end
    xlims!(ax, range1, range2)
    if label && labelpos == :top
        ylims!(ax, 0.875 - height  - (maximum(rows) - 1) * (0.25 + height), 1.375)
    else
        ylims!(ax, 0.975 - height  - (maximum(rows) - 1) * (0.025 + height), 1.05)
    end
    Label(ga[1, 1, Bottom()], "~$(round(range1 / 1e6; digits = 1)) Mb", 
        textsize = 6, halign = :left, valign = :top)
    Label(ga[1, 1, Bottom()], "Chr $(chromosome)", 
        textsize = 6, halign = :center, valign = :top)
    Label(ga[1, 1, Bottom()], "~$(round(range2 / 1e6; digits = 1)) Mb", 
        textsize = 6, halign = :right, valign = :top)
    if label && labelpos == :top
        ax.aspect = AxisAspect(306 / (20 * maximum(rows) * (0.25 + height) / 0.5))
        rowsize!(ga, 1, 18 * maximum(rows) * (0.25 + height) / 0.5)
    elseif label && labelpos != :top
        ax.aspect = AxisAspect(260 / (20 * maximum(rows) * (0.025 + height) / 0.5))
        rowsize!(ga, 1, 18 * maximum(rows) * (0.025 + height) / 0.5)
    else
        ax.aspect = AxisAspect(306 / (20 * maximum(rows) * (0.025 + height) / 0.5))
        rowsize!(ga, 1, 18 * maximum(rows) * (0.025 + height) / 0.5)
    end
    save("$(filename).pdf", f, pt_per_unit = 1)
end