function ellipse(x, a, b)
    if 1 - x^2/a^2 < 0
        return NaN
    end
    abs(b) * sqrt(1 - x^2/a^2)
end
function drawloop!(
    ax::Axis,
    pairedends,
    range1::Real,
    range2::Real;
    height = 100,
    linewidth = 0.25,
    colorarc = "#9658B2",
    colorend = ("#FFBB00", 0.6),
    resolution = 1000 # Plot `resolution` points along x-axis
    )
    midpoints = (sum(pairedends[1])/2, sum(pairedends[2])/2)
    xs = vcat(range(max(midpoints[1], range1), min(midpoints[2], range2);
               step = (range2 - range1) / resolution),
              min(midpoints[2], range2)) # Add last point in case range drops it
    a = (midpoints[2] - midpoints[1]) / 2
    b = height
    ys = ellipse.(xs .- midpoints[1] .- a, a, b)
    lines!(ax, xs, ys; color = colorarc, linewidth = linewidth)
    feet =
    [Polygon([Point2f(pairedend[1], 0),
              Point2f(pairedend[1], -height/20),
              Point2f(pairedend[2], -height/20),
              Point2f(pairedend[2], 0)])
     for pairedend in pairedends]
    poly!(ax, feet; color = colorend)
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
- `linewidth = 0.25`: the line width of the loops' arcs.
- `colorarc = "#9658B2"`: the color of loops' arcs.
- `colorend = "#9658B2"`: the color of loops' ends.
- `resolution = 1000`: plot `resolution` points along x-axis within the given range.
"""
function plotloops!(
        ax::Axis,
        chromosome::AbstractString,
        range1::Real,
        range2::Real,
        loopdf::AbstractDataFrame;
        ymax::Real = 102,
        linewidth::Real = 0.25,
        colorarc = "#9658B2",
        colorend = ("#FFBB00", 0.6),
        resolution = 1000
    )
    loopdf = subset(loopdf,
                    [:chr1, :chr2] =>
                    (chr1, chr2) -> chr1 .== chr2 .== chromosome,
                    [:x1, :y2] =>
                    (start, stop) -> (start .< range2) .&& (stop .> range1)
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
        drawloop!(ax,
                  [(row.x1, row.x2), (row.y1, row.y2)],
                  range1, range2;
                  height = row.b, linewidth = linewidth,
                  colorarc = colorarc, colorend = colorend,
                  resolution = resolution)
    end
    ax.spinewidth = 0.75
    xlims!(ax, range1, range2)
    ylims!(ax, -ymax/20, ymax)
    hidespines!(ax, :t, :r)
    hidedecorations!(ax)
end

plotloops!(ax::Axis, chromosome::AbstractString, bp::Real, loopdf::AbstractDataFrame; window::Real = 1e6, kwargs...) =
    plotloops!(ax, chromosome, bp - window, bp + window, gwas; kwargs...)

function plotloops!(ax::Axis, gene::AbstractString, loopdf::AbstractDataFrame, gencode::DataFrame; window::Real = 1e6, kwargs...)
    chr, start, stop = findgene(gene, gencode)
    plotloops!(ax, chr, start - window, stop + window, loopdf; kwargs...)
end
