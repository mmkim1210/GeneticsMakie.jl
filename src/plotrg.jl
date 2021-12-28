function plotrg!(ax::Axis, r::AbstractMatrix, cols::AbstractVector; 
    p::Union{AbstractMatrix, Nothing} = nothing, circle::Bool = true)
    n = size(r, 1)
    ps = Vector{Polygon}(undef, n^2)
    counter = 1
    for i in 1:n
        for j in 1:n
            ps[counter] = Polygon([
                Point2f(i - 1, 1 - j),
                Point2f(i - 1, -j),
                Point2f(i, -j),
                Point2f(i, 1 - j)
            ])
            counter += 1
        end
    end
    psd = Vector{Polygon}(undef, n)
    for i in 1:n
        psd[i] = Polygon([
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
                csl[counter] = Circle(Point2f(i - 0.5, 0.5 - j), r[i, j] * 0.45)
                csu[counter] = Circle(Point2f(j - 0.5, 0.5 - i), r[j, i] * 0.45)
                colorl[counter] = r[i, j]
                coloru[counter] = r[j, i]
                counter += 1
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
                csl[counter] = Rect(i - 1 + (1 - widthl) / 2, 1 - j - (1 - widthl) / 2, widthl, -widthl)
                csu[counter] = Rect(j - 1 + (1 - widthu) / 2, 1 - i - (1 - widthu) / 2, widthu, -widthu)
                colorl[counter] = r[i, j]
                coloru[counter] = r[j, i]
                counter += 1
            end
        end
    end
    poly!(ax, ps, color = :white, strokewidth = 0.75, strokecolor = "#C3C3C3")
    poly!(ax, psd, color = "#C3C3C3", strokewidth = 0.75, strokecolor = "#C3C3C3")
    poly!(ax, csl, color = colorl, colorrange = (-1, 1), colormap = :RdBu_10)
    poly!(ax, csu, color = coloru, colorrange = (-1, 1), colormap = :RdBu_10)
    ax.xticks = (0.5:(n - 0.5), cols)
    ax.yticks = (-0.5:-1:-(n - 0.5), cols)
    ax.xticklabelsize = 6
    ax.yticklabelsize = 6
    ax.xticklabelpad = 1
    ax.yticklabelpad = 2
    ax.xaxisposition = :top
    ax.yaxisposition = :left
    hidespines!(ax)
    hidedecorations!(ax, ticklabels = false)
    xlims!(ax, 0, n)
    ylims!(ax, -n, 0)
    ax.aspect = DataAspect()
end