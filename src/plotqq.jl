"""
    plotqq!(ax::Axis, gwas::DataFrame; kwargs)
    plotqq!(ax::Axis, P::AbstractVector; kwargs)

Plot QQ plot of `P` values where the expected distribution is the uniform distribution.

Keyword arguments include `xstep::Real` and `ystep::Real` for x and y axes ticks step sizes.
"""
function plotqq!(
    ax::Axis,
    P::AbstractVector; 
    xstep::Real = 1,
    ystep::Real = 2
    )

    n = length(P)
    high = Array{Float64}(undef, n)
    low = Array{Float64}(undef, n)
    for i in 1:n
        high[i] = -log(10, quantile(Beta(i, n + 1 - i), 0.025))
        low[i] = -log(10, quantile(Beta(i, n + 1 - i), 0.975))
    end
    obs = -log.(10, sort(P))
    obs = clamp.(obs, 0.0, -log(10, floatmin(0.0)))
    expect = -log.(10, (collect(1:n) .- 0.5) ./ n)
    xmax = ceil(maximum(expect))
    ymax = ceil(maximum(obs))
    xticks = 0:xstep:xmax
    yticks = (ymax % 2 == 0) ? (0:ystep:ymax) : (0:ystep:(ymax + 1))
    λgc = cquantile(Chisq(1), median(P)) / quantile(Chisq(1), 0.5)

    band!(ax, expect, low, high, color = :gray90)
    scatter!(ax, expect, obs, markersize = 2, color = :black)
    lines!(ax, expect, expect, color = :red2, linewidth = 0.75)
    text!(ax, xmax - 0.5, 1, text = "λgc = $(round(λgc, digits = 2))", 
        fontsize = 6, align = (:right, :bottom))
    xlims!(ax, 0, xmax)
    ylims!(ax, 0, ymax)
    hidedecorations!(ax, label = false, ticklabels = false, ticks = false)
    ax.aspect = AxisAspect(1)
    ax.xlabel = rich("Expected -log", subscript("10"), rich("P", font = :italic))
    ax.ylabel = rich("Observed -log", subscript("10"), rich("P", font = :italic))
    ax.xticklabelsize = 6
    ax.yticklabelsize = 6
    ax.xlabelsize = 8
    ax.ylabelsize = 8
    ax.xticksize = 3
    ax.yticksize = 3
    ax.xtickwidth = 0.75
    ax.ytickwidth = 0.75
    ax.spinewidth = 0.75
    ax.xticks = xticks
    ax.yticks = yticks
end

plotqq!(ax::Axis, gwas::DataFrame; kwargs...) = plotqq!(ax, gwas.P; kwargs...)