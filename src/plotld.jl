"""
    plotld(LD::AbstractMatrix; kwargs)
    plotld!(ax::Axis, LD::AbstractMatrix; kwargs)

Heatmap of symmetric correlation matrix `LD` with the diagonal elements on the x-axis.

# Keyword arguments
```
threshold       threshold below which values are ignored; default 1/9
colormap        colormap of values; default cgrad(:Blues_9, 9, categorical = true)
colorrange      start and end points of colormap; default (0, 1)
strokewidth     width of outline around heatmap boxes; default 0
```
"""
@recipe(PlotLD, LD) do scene
    Attributes(
        threshold   = 1 / 9,
        colormap    = cgrad(:Blues_9, 9, categorical = true),
        colorrange  = (0, 1),
        strokewidth = 0
    )
end

function Makie.plot!(plot::PlotLD{<:Tuple{<:AbstractMatrix}})
    LD = plot[:LD][]
    threshold = plot[:threshold][]
    n = size(LD, 1)
    m = count(>(threshold), LD)
    addx1 = [0, 1, 0, 1]
    addx2 = [1, 1, 0, 1]
    addy  = [0, -1, 0, 1]
    polys = Vector{Vector{Point{2, Float32}}}(undef, Int((m - n) / 2 + n))
    LDvech = Vector{Float32}(undef, Int((m - n) / 2 + n))
    counter = 1
    for i in 1:n
        for j in i:n
            if LD[i, j] <= threshold
                continue
            else
                polys[counter] = [Point2f(
                    abs(i - j) * 5 / n + 10 / n * (i - addx2[k]) + 5 / n * addx1[k], # x coords
                    5 - abs(i - j) * 5 / n + 5 / n * addy[k] # y coords
                    ) for k in 1:4]
                LDvech[counter] = LD[i, j]
                counter += 1
            end
        end
    end
    poly!(plot, polys, color = LDvech, colorrange = plot[:colorrange][],
        colormap = plot[:colormap][], strokewidth = plot[:strokewidth][])
    plot
end
