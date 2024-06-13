"""
    plotrg(r::AbstractMatrix)
    plotrg!(ax::Axis, r::AbstractMatrix)

Correlation plot of matrix `r`.

# Keyword arguments
```
circle          whether to draw cicles instead of rectangles; default true
diagonal        whether to visualize diagonal elements; default false
colormap        colormap of values; default :RdBu_10
colorrange      start and end points of colormap; default (-1, 1)
strokewidth     width of outline around surrounding boxes; default 0.5
```
"""
@recipe(Plotrg, r) do scene
    Attributes(
        circle      = true,
        diagonal    = false,
        colormap    = :RdBu_10,
        colorrange  = (-1, 1),
        strokewidth = 0.5
    )
end

function Makie.plot!(plot::Plotrg{<:Tuple{<:AbstractMatrix}})
    r           = plot[:r][]
    circle      = plot[:circle][]
    diagonal    = plot[:diagonal][]
    colormap    = plot[:colormap][]
    colorrange  = plot[:colorrange][]
    strokewidth = plot[:strokewidth][]
    n = size(r, 1)
    polys = Vector{Polygon}(undef, n^2)
    counter = 1
    for i in 1:n
        for j in 1:n
            polys[counter] = Polygon([
                Point2f(i - 1, 1 - j),
                Point2f(i - 1, -j),
                Point2f(i, -j),
                Point2f(i, 1 - j)
            ])
            counter += 1
        end
    end
    polysd = Vector{Polygon}(undef, n)
    for i in 1:n
        polysd[i] = Polygon([
            Point2f(i - 1, 1 - i),
            Point2f(i - 1, -i),
            Point2f(i, -i),
            Point2f(i, 1 - i)
        ])
    end
    if circle
        csl = Vector{Circle}(undef, binomial(n, 2))
        csu = Vector{Circle}(undef, binomial(n, 2))
        colorl = Vector{Float32}(undef, binomial(n, 2))
        coloru = Vector{Float32}(undef, binomial(n, 2))
        counter = 1
        for i in 1:n
            for j in (i + 1):n
                csl[counter] = Circle(Point2f(j - 0.5, 0.5 - i), r[i, j] * 0.45)
                csu[counter] = Circle(Point2f(i - 0.5, 0.5 - j), r[j, i] * 0.45)
                colorl[counter] = r[i, j]
                coloru[counter] = r[j, i]
                counter += 1
            end
        end
        if diagonal
            csd = Vector{Circle}(undef, n)
            colord = Vector{Float32}(undef, n)
            for i in 1:n
                csd[i] = Circle(Point2f(i - 0.5, 0.5 - i), r[i, i] * 0.45)
                colord[i] = r[i, i]
            end
        end
    else
        csl = Vector{Rect}(undef, binomial(n, 2))
        csu = Vector{Rect}(undef, binomial(n, 2))
        colorl = Vector{Float32}(undef, binomial(n, 2))
        coloru = Vector{Float32}(undef, binomial(n, 2))
        counter = 1
        for i in 1:n
            for j in (i + 1):n
                widthl = r[i, j] * 0.95
                widthu = r[j, i] * 0.95
                csl[counter] = Rect(j - 1 + (1 - widthl) / 2, 1 - i - (1 - widthl) / 2, widthl, -widthl)
                csu[counter] = Rect(i - 1 + (1 - widthu) / 2, 1 - j - (1 - widthu) / 2, widthu, -widthu)
                colorl[counter] = r[i, j]
                coloru[counter] = r[j, i]
                counter += 1
            end
        end
        if diagonal
            csd = Vector{Rect}(undef, n)
            colord = Vector{Float32}(undef, n)
            for i in 1:n
                widthd = r[i, i] * 0.95
                csd[i] = Rect(i - 1 + (1 - widthd) / 2, 1 - i - (1 - widthd) / 2, widthd, -widthd)
                colord[i] = r[i, i]
            end
        end
    end
    poly!(plot, polys, color = :white, strokewidth = strokewidth, strokecolor = "#C3C3C3")
    # poly!(plot, polysd, color = "#C3C3C3", strokewidth = strokewidth, strokecolor = "#C3C3C3")
    poly!(plot, csl, color = colorl, colorrange = colorrange, colormap = colormap)
    poly!(plot, csu, color = coloru, colorrange = colorrange, colormap = colormap)
    if diagonal
        poly!(plot, csd, color = colord, colorrange = colorrange, colormap = colormap)
    end
end