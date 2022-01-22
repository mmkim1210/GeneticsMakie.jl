"""
    coordinategwas(gwas::Vector{DataFrame}; freey::Bool)

Determine shared coordinates for a set of `gwas`.
"""
function coordinategwas(gwas::Vector{DataFrame}; freey::Bool = false)
    df = gwas[1][:, [:CHR, :BP, :P]]
    rename!(df, :P => string("P", 1))    
    for i in 2:length(gwas)
        storage = select(gwas[i], :CHR, :BP, :P)
        rename!(storage, :P => string("P", i))
        df = innerjoin(df, storage; on = [:CHR, :BP])
    end
    for i in 1:length(gwas)
        df[!, string("log10P", i)] = -log.(10, df[!, string("P", i)])
    end
    if !freey
        ymax = 0
        for i in 1:length(gwas)
            storage = ceil(maximum(df[!, string("log10P", i)])) + 2.5
            ymax = max(ymax, storage)
        end
        ymaxs = fill(ymax, length(gwas))
    else
        ymaxs = Vector{Float64}(undef, length(gwas))
        for i in 1:length(gwas)
            ymaxs[i] = ceil(maximum(df[:, string("log10P", i)])) + 2.5
        end
    end
    dforder = DataFrame("CHR" => vcat(string.(1:22), ["X", "Y"]), "rank" => 1:24)
    storage = combine(groupby(df, :CHR), :BP => maximum => :maxpos)
    storage = leftjoin(storage, dforder; on = :CHR)
    sort!(storage, order(:rank))
    storage.add = cumsum(storage.maxpos) - storage.maxpos
    df = leftjoin(df, storage; on = :CHR)
    df.x = df.BP + df.add
    xmax = maximum(df.x)
    ticks = combine(groupby(df, :CHR), :x => middle => :center)
    return df, ymaxs, xmax, ticks
end

"""
    plotgwas!(ax::Axis, df::DataFrame, i::Int, ymax::Real, xmax::Real, ticks::DataFrame; ystep::Real, sigline::Bool, sigcolor::Bool)

Plot gwas results with coordinates from `coordinategwas`, namely `df` for phenotype `i` with
x and y limits, `xmax` and `ymax`, respectively. The position of x ticks are determined by `ticks`.
"""
function plotgwas!(
    ax::Axis,
    df::DataFrame,
    i::Int,
    ymax::Real,
    xmax::Real,
    ticks::DataFrame; 
    ystep::Real = 2,
    sigline::Bool = false,
    sigcolor::Bool = true
)

    indeven = findall(in(vcat(string.(1:2:22), "X")), df.CHR)
    indodd = findall(in(vcat(string.(2:2:23), "Y")), df.CHR)
    scatter!(ax, view(df, indeven, :x), view(df, indeven, string("log10P", i)), 
        markersize = 1.5, color = "#0D0D66")
    scatter!(ax, view(df, indodd, :x), view(df, indodd, string("log10P", i)), 
        markersize = 1.5, color = "#7592C8")
    if sigcolor
        ind = df[!, string("log10P", i)] .> -log(10, 5e-8)
        dfsig = view(df, ind, :)
        scatter!(ax, dfsig.x, dfsig[!, string("log10P", i)], 
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
    xticklabels = unique(ticks.CHR)
    ax.xlabel = "Chromosome"
    ax.ylabel = "-log[p]"
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