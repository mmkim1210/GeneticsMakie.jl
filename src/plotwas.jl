function plotwas!(
    ax::Axis,
    twas::DataFrame,
    gencode::DataFrame; 
    ymax::Real = 0,
    sigline::Bool = false,
    sigcolor::Bool = true,
    zscore::Bool = false,
    gene_id::Bool = false
)
    gencodeₛ = filter(x -> x.feature == "gene", gencode)
    gencodeₛ.x = (gencodeₛ.start + gencodeₛ.end) / 2
    if gene_id
        select!(gencodeₛ, :gene_id, :x, :CHR)
        df = innerjoin(twas, gencodeₛ; on = :gene_id)
    else
        select!(gencodeₛ, :gene_name, :x, :CHR)
        df = innerjoin(twas, gencodeₛ; on = :gene_name)
    end
    df = select(gwas, :CHR, :BP, :P)
    df.P = -log.(10, df.P)
    if ymax == 0
        ymax = maximum(df.P) / 4 * 5
        ymax <= 10 ? ymax = 10 : nothing
    end
    storage = DataFrame(CHR = vcat(string.(1:22), ["X", "Y"]))
    storage.maxpos = [GRCh37_totlength[chr] for chr in storage.CHR]
    storage.add = cumsum(storage.maxpos) - storage.maxpos
    df = leftjoin(df, storage; on = :CHR)
    df.x = df.BP + df.add
    xmax = sum(unique(df.maxpos))
    indeven = findall(in(vcat(string.(1:2:22), "X")), df.CHR)
    indodd = findall(in(vcat(string.(2:2:23), "Y")), df.CHR)
    scatter!(ax, view(df, indeven, :x), view(df, indeven, :P), markersize = 1.5, color = "#0D0D66")
    scatter!(ax, view(df, indodd, :x), view(df, indodd, :P), markersize = 1.5, color = "#7592C8")
    if sigcolor
        ind = df.P .> -log(10, 5e-8)
        dfsig = view(df, ind, :)
        scatter!(ax, dfsig.x, dfsig.P,  markersize = 1.5, color = "#4DB069")
    end
    if sigline
        hlines!(ax, -log(10, 5e-8); xmin = 0.0, xmax = xmax, linewidth = 0.75, color = :red2)
    end
    xlims!(ax, 0, xmax)
    ylims!(ax, 0, ymax)
    hideydecorations!(ax, label = false, ticklabels = false, ticks = false)
    hidexdecorations!(ax, label = false, ticklabels = false)
    ax.xlabel = "Chromosome"
    ax.ylabel = "-log[p]"
    ax.xticks = ((cumsum(storage.maxpos) + storage.add) / 2, storage.CHR)
    ax.yticks = setticks(ymax)
    ax.xticklabelsize = 6
    ax.yticklabelsize = 6
    ax.xlabelsize = 8
    ax.ylabelsize = 8
    ax.xticksize = 3
    ax.yticksize = 3
    ax.xtickwidth = 0.75
    ax.ytickwidth = 0.75
    ax.spinewidth = 0.75
    ax.xticklabelpad = 0
    ax.yticklabelpad = 2.5
end

function plotwas!(
    ax::Axis,
    chromosome::AbstractString,
    range1::Real,
    range2::Real,
    twas::DataFrame,
    gencode::DataFrame
)

end

plotwas!(ax::Axis, chromosome::AbstractString, bp::Real, twas::DataFrame, gencode::DataFrame; window::Real = 1e6) =
    plotwas!(ax, chromosome, bp - window, bp + window, twas, gencode)

function plotwas!(ax::Axis, gene::AbstractString, twas::DataFrame, gwas::DataFrame; window::Real = 1e6)
    chr, start, stop = findgene(gene, gencode)
    plotwas!(ax, chr, start - window, stop + window, twas, gencode)
end


# g <- ggplot() + geom_point(data=df, aes(x=pos, y=TWAS.Z, col=(CHR %% 2 == 0))) + 
#   geom_point(data=df[df$FOCUS==TRUE,], aes(x=pos, y=TWAS.Z, col="3")) + 
#   scale_color_manual(values=c("FALSE"="gray60","TRUE"="gray70","3"="#E41A1C")) +
#   scale_x_continuous("Chromosome", breaks=floor(0.5*(ticks+ticks2)), labels=1:22) +
#   geom_text_repel(data=df[df$FOCUS==TRUE,], aes(x=pos, y=TWAS.Z, label=ID), 
#                   color="black", size=3, segment.size=.1) + 
#   theme_bw() + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank(), 
#                      panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), legend.position="none") + 
#   ylim(-10,10) + ylab("TWAS Z-score") + geom_hline(yintercept=4.63, linetype="dashed", size=0.375, color="#377EB8") + 
#   geom_hline(yintercept=-4.63, linetype="dashed", size=0.375, color="#377EB8") + 
#   geom_vline(xintercept=ticks2-1, linetype="dotted", size=0.375, color="black")

function coordinateisforms(
    gene::AbstractString, 
    gencode::DataFrame,
    orderby::Union{Nothing, AbstractVector{<:AbstractString}},
    height::Real,
    text::Union{Bool, Symbol},
    shrink::Real
    )

    df = filter(x -> x.gene_name == gene, gencode)
    nrow(df) == 0 ? error("Cannot find $(gene) in the annotation.") : nothing
    dfi = view(df, df.feature .== "transcript", :)
    dfe = view(df, df.feature .== "exon", :)
    chromosome = df.seqnames[1]
    isoforms = unique(dfi.transcript_id)
    n = length(isoforms)
    ps = Vector{Vector{Polygon}}(undef, n) 
    bs = Matrix{Float64}(undef, n, 2)
    rows = collect(1:n)
    if shrink != 1
        for (exon1, exon2) in zip(storage.start[2:end], storage.end[1:end])
            (exon1 - exon2) * shrink
        end
    end
    if !isnothing(orderby)
        prior = []
        for k in 1:length(orderby)
            if any(isoforms .== orderby[k])
                push!(prior, orderby[k])
            end
        end
        if length(prior) > 0
            dforder = DataFrame(transcript_id = [prior; filter(!in(orderby), isoforms)], rank = 1:n)
            dfi = leftjoin(dfi, dforder; on = :transcript_id)
            sort!(dfi, [order(:rank)])
            isoforms = unique(dfi.transcript_id)
        end
    end
    for j in eachindex(isoforms)
        ranges = view(dfe, findall(isequal(isoforms[j]), dfe.transcript_id), [:start, :end])
        m = size(ranges, 1)
        p = Vector{Polygon}(undef, m)
        if text == true || text == :top || text == :t || text == :bottom || text == :b
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
    plotisoforms!(ax::Axis, gene::AbstractString, gencode::DataFrame; kwargs)

Plot each isoform of a given `gene` on a separate row.

# Arguments
- `orderby::Union{Nothing, AbstractVector{<:AbstractString}} = nothing`: the order of isoforms.
- `highlight::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing`: isoforms to be highlighted and their colors.
- `height::Real = 0.25`: the height of exons.
- `isoformcolor = :royalblue`: the color of isoforms.
- `textcolor = :black`: the color of isoform labels.
- `text::Union{Bool, Symbol} = :top`: the position of isoform labels. 
"""
function plotisoforms!(
    ax::Axis,
    gene::AbstractString,
    gencode::DataFrame;
    orderby::Union{Nothing, AbstractVector{<:AbstractString}} = nothing,
    highlight::Union{Nothing, Tuple{AbstractVector, AbstractVector}} = nothing, 
    height::Real = 0.25,
    isoformcolor = :royalblue,
    textcolor = :black,
    text::Union{Bool, Symbol} = :top,
    shrink::Real = 1
    )

    isoforms, ps, bs, rows, chromosome = coordinateisforms(gene, gencode, orderby, height, text, shrink)
    isnothing(highlight) ? highlight = ([nothing], [nothing]) : nothing
    if text == true || text == :top || text == :t
        for j in 1:size(ps, 1)
            ind = findfirst(isequal(isoforms[j]), highlight[1])
            if isnothing(ind)
                poly!(ax, ps[j], color = isoformcolor, strokewidth = 0)
                lines!(ax, [bs[j, 1], bs[j, 2]], 
                    [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)],
                    color = isoformcolor, linewidth = 0.5)
                text!(ax, "$(isoforms[j])", 
                    position = ((bs[j, 1] + bs[j, 2]) / 2, 1 - (rows[j] - 1) * (0.25 + height)), 
                    align = (:center, :bottom), fontsize = 6, color = textcolor)
            else
                poly!(ax, ps[j], color = highlight[2][ind], strokewidth = 0)
                lines!(ax, [bs[j, 1], bs[j, 2]], 
                    [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)],
                    color = highlight[2][ind], linewidth = 0.5)
                text!(ax, "$(isoforms[j])", 
                    position = ((bs[j, 1] + bs[j, 2]) / 2, 1 - (rows[j] - 1) * (0.25 + height)), 
                    align = (:center, :bottom), fontsize = 6, color = textcolor)
            end
        end
    elseif text == :bottom || text == :b
        for j in 1:size(ps, 1)
            ind = findfirst(isequal(isoforms[j]), highlight[1])
            if isnothing(ind)
                poly!(ax, ps[j], color = isoformcolor, strokewidth = 0)
                lines!(ax, [bs[j, 1], bs[j, 2]], 
                    [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)],
                    color = isoformcolor, linewidth = 0.5)
                text!(ax, "$(isoforms[j])", 
                    position = ((bs[j, 1] + bs[j, 2]) / 2, 1 - height - (rows[j] - 1) * (0.25 + height)), 
                    align = (:center, :top), fontsize = 6, color = textcolor)
            else
                poly!(ax, ps[j], color = highlight[2][ind], strokewidth = 0)
                lines!(ax, [bs[j, 1], bs[j, 2]], 
                    [1 - height / 2 - (rows[j] - 1) * (0.25 + height), 1 - height / 2 - (rows[j] - 1) * (0.25 + height)],
                    color = highlight[2][ind], linewidth = 0.5)
                text!(ax, "$(isoforms[j])", 
                    position = ((bs[j, 1] + bs[j, 2]) / 2, 1 - height - (rows[j] - 1) * (0.25 + height)), 
                    align = (:center, :top), fontsize = 6, color = textcolor)
            end
        end
    else
        for j in 1:size(ps, 1)
            ind = findfirst(isequal(isoforms[j]), highlight[1])
            if isnothing(ind)
                poly!(ax, ps[j], color = isoformcolor, strokewidth = 0)
                lines!(ax, [bs[j, 1], bs[j, 2]], 
                    [1 - height / 2 - (rows[j] - 1) * (0.025 + height), 1 - height / 2 - (rows[j] - 1) * (0.025 + height)],
                    color = isoformcolor, linewidth = 0.5)
            else
                poly!(ax, ps[j], color = highlight[2][ind], strokewidth = 0)
                lines!(ax, [bs[j, 1], bs[j, 2]], 
                    [1 - height / 2 - (rows[j] - 1) * (0.025 + height), 1 - height / 2 - (rows[j] - 1) * (0.025 + height)],
                    color = highlight[2][ind], linewidth = 0.5)
            end
        end
    end
    range1 = minimum(bs[:, 1])
    range2 = maximum(bs[:, 2])
    diff = range2 - range1
    prop = 23000 * diff / 2.2e6
    if text == true || text == :top || text == :t || text == :bottom || text == :b
        storage = bs[:, 1] + bs[:, 2]
        ind = argmin(storage)
        low = storage[ind] / 2 - prop * 12
        low < range1 ? range1 = low : range1 -= prop * 4.5
        ind = argmax(storage)
        high = storage[ind] / 2 + prop * 12
        high > range2 ? range2 = high : range2 += prop * 4.5
    else 
        range1 = range1 - diff / 50
        range2 = range2 + diff / 50
    end
    ax.spinewidth = 0.75
    hidexdecorations!(ax)
    if text == :left || text == :l || text == :right || text == :r
        ax.yticks = ([1 - height / 2 - (rows[j] - 1) * (0.025 + height) for j in 1:length(isoforms)], isoforms)
        hideydecorations!(ax, ticklabels = false)
        ax.yticklabelsize = 4
        ax.yticklabelpad = 3
        (text == :right) || (text == :r) ? ax.yaxisposition = :right : ax.yaxisposition = :left
        rs = 18 * maximum(rows) * (0.025 + height) / 0.5
    elseif text == true || text == :top || text == :t || text == :bottom || text == :b
        hideydecorations!(ax)
        rs = 18 * maximum(rows) * (0.25 + height) / 0.5
    else
        hideydecorations!(ax)
        rs = 18 * maximum(rows) * (0.025 + height) / 0.5
    end
    xlims!(ax, range1, range2)
    if text == true || text == :top || text == :t
        ylims!(ax, 0.875 - height  - (maximum(rows) - 1) * (0.25 + height), 1.375)
    elseif text == :bottom || text == :b
        ylims!(ax, 0.65 - height  - (maximum(rows) - 1) * (0.25 + height), 1.05)
    else
        ylims!(ax, 0.975 - height  - (maximum(rows) - 1) * (0.025 + height), 1.05)
    end
    return rs, chromosome, range1, range2
end