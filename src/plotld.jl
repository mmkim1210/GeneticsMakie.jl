"""
    plotld(LD::AbstractMatrix)

Visualize a correlation matrix `LD` with the diagonal elements on the x-axis.
"""
function plotld(LD::AbstractMatrix;
    filename::AbstractString = "LD",
    colorbar = true,
    xlabel::Union{Nothing, Tuple{AbstractString, Real, Real}} = nothing)

    CairoMakie.activate!(type = "png")
    set_theme!(font = "Arial")

    n = size(LD, 1)
    addx1 = [0, 1, 0, 1]
    addx2 = [1, 1, 0, 1]
    addy = [0, -1, 0, 1]

    ps = Vector{Polygon}(undef, binomial(n, 2) + n) # collect coordinates for squares
    LDvech = Vector{Float64}(undef, binomial(n, 2) + n) # collect off-diagonal elements
    counter = 1
    for i in 1:n
        for j in i:n
            ps[counter] = Polygon([Point2f(
                abs(i - j) * 5 / n + 10 / n * (i - addx2[k]) + 5 / n * addx1[k], # x coordinates
                5 - abs(i - j) * 5 / n + 5 / n * addy[k] # y coordinates
                ) for k in 1:4])
            LDvech[counter] = LD[i, j]
            counter += 1
        end
    end
    
    f = Figure(resolution = (306, 200))
    ga = f[1, 1] = GridLayout()
    ax = Axis(ga[1, 1])
    poly!(ax, ps, color = LDvech, colorrange = (0, 1),
        colormap = cgrad(ColorSchemes.Blues_9, 9, categorical = true), strokewidth = 0)
    ax.aspect = DataAspect()
    xlims!(0, 10)
    ylims!(0, 5)
    hidedecorations!(ax)
    if !isnothing(xlabel)
        hidespines!(ax, :r, :l, :b)
        ax.spinewidth = 0.75
        Label(ga[1, 1, Top()], "~$(round(xlabel[2] / 1e6; digits = 1)) Mb",
            textsize = 6, halign = :left, valign = :top)
        Label(ga[1, 1, Top()], "Chr $(xlabel[1])",
            textsize = 6, halign = :center, valign = :top)
        Label(ga[1, 1, Top()], "~$(round(xlabel[3] / 1e6; digits = 1)) Mb",
            textsize = 6, halign = :right, valign = :top)
    else
        hidespines!(ax)
    end
    rowsize!(ga, 1, 137)
    if colorbar
        Colorbar(ga[2, 1], limits = (0, 1), 
            colormap = cgrad(ColorSchemes.Blues_9, 9, categorical = true),
            label = "LD", labelsize = 6, vertical = false, flipaxis = false,
            width = 40, tellwidth = false, ticklabelsize = 6,
            spinewidth = 0.5, tickwidth = 0.5, height = 5, ticksize = 3)
    end
    rowgap!(ga, 5)
    save("$(filename).png", f, px_per_unit = 4)
end