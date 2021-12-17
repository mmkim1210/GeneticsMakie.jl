"""
    coordinategwas(gwas::Vector{DataFrame}; ymax, chr, bp, p, freey)

Determine shared coordinates for a set of `gwas`.
"""
function coordinategwas(gwas::Vector{DataFrame};
    ymax::Real = 0,
    chr::AbstractString = "CHR",
    bp::AbstractString = "BP",
    p::AbstractString = "P", 
    freey::Bool = false)

    for i in 1:length(gwas)
        (all(gwas[i][!, p] .>= 0.0) && all(gwas[i][!, p] .<= 1.0)) || 
            @error("P values are not between 0 and 1 for phenotype $(i).")
    end

    df = gwas[1][:, [chr, bp, p]]
    rename!(df, p => string(p, 1))    
    for i in 2:length(gwas)
        df = @chain gwas[i] begin
            select(chr, bp, p)
            rename(p => string(p, i))
            innerjoin(df, _; on = [chr, bp])
        end
    end
    for i in 1:length(gwas)
        df[!, string("log10", p, i)] = -log.(10, df[!, string(p, i)])
        df[!, string("log10", p, i)] = clamp.(df[!, string("log10", p, i)], 0.0, -log(10, floatmin(0.0)))
    end
    if ymax == 0
        if !freey
            for i in 1:length(gwas)
                storage = ceil(maximum(df[:, string("log10", p, i)])) + 2.5
                ymax = max(ymax, storage)
            end
            ymaxs = fill(ymax, length(gwas))
        else
            ymaxs = Vector{Float64}(undef, length(gwas))
            for i in 1:length(gwas)
                ymaxs[i] = ceil(maximum(df[:, string("log10", p, i)])) + 2.5
            end
        end
    end
    dforder = DataFrame(chr => vcat(string.(1:22), ["X", "Y"]), "rank" => 1:24)
    df = @chain df begin
        groupby(chr)
        combine(bp => maximum => :maxpos)
        leftjoin(_, dforder; on = chr)
        sort(order(:rank))
        @transform(:add = cumsum(:maxpos) - :maxpos)
        leftjoin(df, _; on = chr)
    end
    df[!, :x] = df[!, bp] + df[!, :add]
    xmax = maximum(df.x)
    ticks = @chain df begin
        groupby(chr)
        @combine(:center = middle(:x))
    end
    return df, ymaxs, xmax, ticks
end

"""
    plotgwas!(ax::Axis, df::DataFrame, i::Int, ymax::Real, xmax::Real, ticks::DataFrame; xlabel, ylabel, ystep, chr, p, sigline, sigcolor)

Plot gwas results with coordinates from `coordinategwas`, namely `df` for phenotype `i` with
x and y limits, `xmax` and `ymax`, respectively. The position of x ticks are
determined by `ticks`.
"""
function plotgwas!(ax::Axis,
    df::DataFrame,
    i::Int,
    ymax::Real,
    xmax::Real,
    ticks::DataFrame; 
    xlabel::AbstractString = "Chromosome",
    ylabel::AbstractString = "-log[p]",
    ystep::Real = 2,
    chr::AbstractString = "CHR",
    p::AbstractString = "P",
    sigline::Bool = false,
    sigcolor::Bool = true)

    indeven = findall(in(vcat(string.(1:2:22), "X")), df[!, chr])
    indodd = findall(in(vcat(string.(2:2:23), "Y")), df[!, chr])
    scatter!(ax, df[indeven, :x], df[indeven, string("log10", p, i)], 
        markersize = 1.5, color = "#0D0D66")
    scatter!(ax, df[indodd, :x], df[indodd, string("log10", p, i)], 
        markersize = 1.5, color = "#7592C8")
    if sigcolor
        dfsig = df[df[!, string("log10", p, i)] .> -log(10, 5e-8), :]
        scatter!(ax, dfsig[!, :x], dfsig[!, string("log10", p, i)], 
            markersize = 1.5, color = "#4DB069")
    end
    if sigline
        hlines!(ax, -log(10, 5e-8); xmin = 0.0, xmax = xmax, 
            linewidth = 0.75, color = :red2)
    end
    xlims!(ax, 0, xmax)
    ylims!(ax, 0, ymax)
    hideydecorations!(ax, label = false, ticklabels = false, ticks = false)
    hidexdecorations!(ax, label = false, ticklabels = false)
    xticks = ticks.center
    xticklabels = unique(ticks[!, chr])
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.xticks = (xticks, xticklabels)
    ax.yticks = 0:ystep:ymax
    ax.xticklabelsize = 6
    ax.yticklabelsize = 6
    ax.xlabelsize = 8
    ax.ylabelsize = 8
    ax.xticksize = 3
    ax.yticksize = 3
    ax.xtickwidth = 0.75
    ax.ytickwidth = 0.75
    ax.spinewidth = 0.75
    ax.xticklabelpad = 0
    ax.yticklabelpad = 2.5
end