"""
    plotgwas!(ax::Axis, gwas::DataFrame; ymax::Real, sigline::Bool, sigcolor::Bool)

Plot `gwas` results as a Manhattan plot.
"""
function plotgwas!(
    ax::Axis,
    gwas::DataFrame; 
    ymax::Real = 0,
    sigline::Bool = false,
    sigcolor::Bool = true
    )

    df = select(gwas, :CHR, :BP, :P)
    df.P = -log.(10, df.P)
    if ymax == 0
        ymax = maximum(df.P) / 4 * 5
        ymax <= 10 ? ymax = 10 : nothing
    end
    storage = DataFrame(CHR = vcat(string.(1:22), ["X", "Y"]))
    storage.maxpos = [GRCh37_totlength[chr] for chr in storage.CHR]
    storage.add = cumsum(storage.maxpos) - storage.maxpos
    df = leftjoin(df, storage; on = :CHR)
    df.x = df.BP + df.add
    xmax = sum(unique(df.maxpos))
    indeven = findall(in(vcat(string.(1:2:22), "X")), df.CHR)
    indodd = findall(in(vcat(string.(2:2:23), "Y")), df.CHR)
    scatter!(ax, view(df, indeven, :x), view(df, indeven, :P), markersize = 1.5, color = "#0D0D66")
    scatter!(ax, view(df, indodd, :x), view(df, indodd, :P), markersize = 1.5, color = "#7592C8")
    if sigcolor
        ind = df.P .> -log(10, 5e-8)
        dfsig = view(df, ind, :)
        scatter!(ax, dfsig.x, dfsig.P,  markersize = 1.5, color = "#4DB069")
    end
    if sigline
        hlines!(ax, -log(10, 5e-8); xmin = 0.0, xmax = xmax, linewidth = 0.75, color = :red2)
    end
    xlims!(ax, 0, xmax)
    ylims!(ax, 0, ymax)
    hideydecorations!(ax, label = false, ticklabels = false, ticks = false)
    hidexdecorations!(ax, label = false, ticklabels = false)
    ax.xlabel = "Chromosome"
    ax.ylabel = "-log[p]"
    ax.xticks = ((cumsum(storage.maxpos) + storage.add) / 2, storage.CHR)
    ax.yticks = setticks(ymax)
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