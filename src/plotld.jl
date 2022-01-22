"""
    plotld!(ax::Axis, LD::AbstractMatrix; color::AbstractString)

Visualize a symmetric correlation matrix `LD` with the diagonal elements on the x axis.
`color` can be either "red", "green", "blue", or "black".
"""
function plotld!(ax::Axis, LD::AbstractMatrix; color::AbstractString = "blue")
    n = size(LD, 1)
    m = count(>(1 / 9), LD)
    addx1 = [0, 1, 0, 1]
    addx2 = [1, 1, 0, 1]
    addy = [0, -1, 0, 1]
    ps = Vector{Polygon}(undef, Int((m - n) / 2 + n))
    LDvech = Vector{Float32}(undef, Int((m - n) / 2 + n))
    counter = 1
    for i in 1:n
        for j in i:n
            if LD[i, j] <= 1 / 9
                continue
            else
                ps[counter] = Polygon([Point2f(
                    abs(i - j) * 5 / n + 10 / n * (i - addx2[k]) + 5 / n * addx1[k], # x coords
                    5 - abs(i - j) * 5 / n + 5 / n * addy[k] # y coords
                    ) for k in 1:4])
                LDvech[counter] = LD[i, j]
                counter += 1
            end
        end
    end
    if color == "black"
        poly!(ax, ps, color = LDvech, colorrange = (0, 1),
            colormap = cgrad(:Greys_9, 9, categorical = true), strokewidth = 0)
    elseif color == "green"
        poly!(ax, ps, color = LDvech, colorrange = (0, 1),
            colormap = cgrad(:Greens_9, 9, categorical = true), strokewidth = 0)
    elseif color == "red"
        poly!(ax, ps, color = LDvech, colorrange = (0, 1),
            colormap = cgrad(:Reds_9, 9, categorical = true), strokewidth = 0)
    else
        poly!(ax, ps, color = LDvech, colorrange = (0, 1),
            colormap = cgrad(:Blues_9, 9, categorical = true), strokewidth = 0)
    end
    ax.aspect = DataAspect()
    ax.spinewidth = 0.75
    xlims!(0, 10)
    ylims!(0, 5)
    hidedecorations!(ax)
    hidespines!(ax, :r, :l, :b)
end