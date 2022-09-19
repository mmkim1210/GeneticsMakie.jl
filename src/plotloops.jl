function ellipse(x, a, b)
    if 1 - x^2/a^2 < 0
        return 0
    end
    abs(b) * sqrt(1 - x^2/a^2)
end
function ellipseband!(
    ax::Axis,
    outerranges,
    innerranges,
    range1::Real,
    range2::Real;
    height = 100,
    color = "#9658B2",
    step = nothing, # Plot a point for each <step> along x-axis
    length = 1000, # Overrides step, sets how many plotted points
    outline = false
    )
    if isnothing(length)
        xs = range(max(outerranges[1], range1), min(outerranges[2], range2);
                   step = step)
    else
        xs = range(max(outerranges[1], range1), min(outerranges[2], range2);
                   length = length)
    end
    outera = (outerranges[2] - outerranges[1]) / 2
    outerb = height
    yhighs = ellipse.(xs .- outerranges[1] .- outera, outera, outerb)
    innera = (innerranges[2] - innerranges[1]) / 2
    innerb = innera / outera * height
    ylows = ellipse.(xs .- innerranges[1] .- innera, innera, innerb)
    band!(ax, xs, ylows, yhighs; color = color)
    if outline
        lines!(ax, xs, yhighs; color = color, linewidth = 0.1)
        lines!(ax, xs, ylows; color = color, linewidth = 0.1)
    end
end

"""
    plotloops!(ax::Axis, chromosome::AbstractString, range1::Real, range2::Real, loopdf::DataFrame; kwargs)
    plotloops!(ax::Axis, chromosome::AbstractString, bp::Real, loopdf::DataFrame; kwargs)
    plotloops!(ax::Axis, gene::AbstractString, loopdf::DataFrame, gencode::DataFrame; kwargs)

Plot loops present in `loopdf` within a given `chromosome` and genomic range between `range1` 
and `range2`.

Alternatively, plot within a given `chromosome` and a certain `window` around a 
genomic coordinate `bp` or plot within a certain `window` around `gene`.

# Arguments
- `ymax::Real = 102`: the maximum value for y axis. 
- `outline::Bool = false`: draw outline around loops (useful for 2.5+ Mbase windows)
- `color = "#9658B2"`: the color of loops.
- `step = nothing`: plot a point for each `step` base pairs. `length` must be set to nothing.
- `length = 1000`: plot `length` number of points for each loop. overrides `step`.
"""
function plotloops!(
        ax::Axis,
        chromosome::AbstractString,
        range1::Real,
        range2::Real,
        loopdf::AbstractDataFrame;
        ymax::Real = 102,
        outline::Bool = false,
        color = "#9658B2",
        step = nothing,
        length = 1000
    )
    loopdf = subset(loopdf,
                    [:chr1, :chr2] =>
                    (chr1, chr2) -> chr1 .== chr2 .== chromosome,
                    [:x1, :x2, :y1, :y2] =>
                    ByRow((coords...) -> any(range1 .< coords .< range2))
                   )
    transform!(loopdf,
               [:x1, :y2] =>
               ((coord1, coord2) -> coord2 .- coord1) =>
               :dist)
    sort!(loopdf, :dist)
    if nrow(loopdf) == 1
        loopdf.b = [100]
    else
        loopdf.b = range(10, 100, length = nrow(loopdf))
    end
    for row in eachrow(loopdf)
        ellipseband!(ax,
                     (row.x1, row.y2), (row.x2, row.y1),
                     range1, range2; height = row.b, color = color,
                     step = step, length = length, outline = outline)
    end
    ax.spinewidth = 0.75
    xlims!(ax, range1, range2)
    ylims!(ax, 0, ymax)
    hidespines!(ax, :t, :r)
    hidedecorations!(ax)
end

plotloops!(ax::Axis, chromosome::AbstractString, bp::Real, loopdf::AbstractDataFrame; window::Real = 1e6, kwargs...) =
    plotloops!(ax, chromosome, bp - window, bp + window, gwas; kwargs...)

function plotloops!(ax::Axis, gene::AbstractString, loopdf::AbstractDataFrame, gencode::DataFrame; window::Real = 1e6, kwargs...)
    chr, start, stop = findgene(gene, gencode)
    plotloops!(ax, chr, start - window, stop + window, loopdf; kwargs...)
end
