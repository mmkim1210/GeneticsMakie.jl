"""
    plotqq!(ax::Axis, p::AbstractVector; title, xlabel, ylabel, xstep, ystep)

Plot QQ plot of `p` values where the expected distribution is the uniform distribution. 
"""
function plotqq!(ax::Axis,
    p::AbstractVector; 
    title::AbstractString = "",
    xlabel::AbstractString = "Expected -log[p]",
    ylabel::AbstractString = "Observed -log[p]",
    xstep::Real = 1,
    ystep::Real = 2)

    (all(p .>= 0.0) && all(p .<= 1.0)) || @error("P values are not between 0 and 1.")

    n = length(p)
    high = Array{Float64}(undef, n)
    low = Array{Float64}(undef, n)
    for i in 1:n
        high[i] = -log(10, quantile(Beta(i, n + 1 - i), 0.025))
        low[i] = -log(10, quantile(Beta(i, n + 1 - i), 0.975))
    end
    obs = -log.(10, sort(p))
    obs = clamp.(obs, 0.0, -log(10, floatmin(0.0)))
    expect = -log.(10, (collect(1:n) .- 0.5) ./ n)
    xmax = ceil(maximum(expect))
    ymax = ceil(maximum(obs))
    xticks = 0:xstep:xmax
    yticks = (ymax % 2 == 0) ? (0:ystep:ymax) : (0:ystep:(ymax + 1))
    λgc = cquantile(Chisq(1), median(p)) / quantile(Chisq(1), 0.5)

    f = Figure()
    band!(ax, expect, low, high, color = :gray90)
    scatter!(ax, expect, obs, markersize = 2, color = :black)
    lines!(ax, expect, expect, color = :red2, linewidth = 0.75)
    text!(ax, "λgc = $(round(λgc, digits = 2))", position = (xmax - 0.5, 1), 
        textsize = 6, align = (:right, :bottom))
    xlims!(ax, 0, xmax)
    ylims!(ax, 0, ymax)
    hidedecorations!(ax, label = false, ticklabels = false, ticks = false)
    ax.aspect = AxisAspect(1)
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.title = title
    ax.xticklabelsize = 6
    ax.yticklabelsize = 6
    ax.xlabelsize = 8
    ax.ylabelsize = 8
    ax.xticksize = 3
    ax.yticksize = 3
    ax.xtickwidth = 0.75
    ax.ytickwidth = 0.75
    ax.spinewidth = 0.75
    ax.titlesize = 8
    ax.xticks = xticks
    ax.yticks = yticks
end

function plotqq!(ax::Axis, df::DataFrame; p::AbstractString = "P", kwargs...)
    if !(in(p, names(df)))
        throw(ArgumentError("column name :" * p * " not found in the dataframe."))
    end
    plotqq!(ax, df[!, p]; kwargs...)
end