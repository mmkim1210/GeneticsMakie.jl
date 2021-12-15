"""
    plotld!(LD::AbstractMatrix)

Visualize a correlation matrix `LD` with the diagonal elements on the x-axis.
"""
function plotld!(ax::Axis, LD::AbstractMatrix)
    n = size(LD, 1)
    addx1 = [0, 1, 0, 1]
    addx2 = [1, 1, 0, 1]
    addy = [0, -1, 0, 1]
    ps = Vector{Polygon}(undef, binomial(n, 2) + n)
    LDvech = Vector{Float32}(undef, binomial(n, 2) + n)
    counter = 1
    for i in 1:n
        for j in i:n
            ps[counter] = Polygon([Point2f(
                abs(i - j) * 5 / n + 10 / n * (i - addx2[k]) + 5 / n * addx1[k], # x coords
                5 - abs(i - j) * 5 / n + 5 / n * addy[k] # y coords
                ) for k in 1:4])
            LDvech[counter] = LD[i, j]
            counter += 1
        end
    end
    poly!(ax, ps, color = LDvech, colorrange = (0, 1),
        colormap = cgrad(ColorSchemes.Blues_9, 9, categorical = true), strokewidth = 0)
    ax.aspect = DataAspect()
    ax.spinewidth = 0.75
    xlims!(0, 10)
    ylims!(0, 5)
    hidedecorations!(ax)
    hidespines!(ax, :r, :l, :b)
end