"""
    plotgwas!(ax::Axis, gwas::DataFrame)

Plot `gwas` results as a Manhattan plot.

# Keyword arguments
```
ymax            maximum value for y-axis
p               genome-wide significance threshold; default 5e-8
linecolor       color of genome-wide significance line, which can be turned off by setting to nothing; default :red2
scattercolor    color of genome-wide significant variants, which can be turned off by setting to nothing; default "#4DB069"
chromcolors     colors of even and odd chromosomes; default ["#0D0D66", "#7592C8"]
build::Int      human genome build; default 37 
```

"""
function plotgwas!(
    ax::Axis,
    gwas::DataFrame; 
    ymax::Real = 0,
    p::Real = 5e-8,
    linecolor = :red2,
    scattercolor = "#4DB069",
    chromcolors = ["#0D0D66", "#7592C8"],
    build::Int = 37
    )

    df = select(gwas, :CHR, :BP, :P)
    df.P = -log.(10, df.P)
    if ymax == 0
        ymax = maximum(df.P) / 4 * 5
        ymax <= 10 ? ymax = 10 : nothing
    end
    storage = DataFrame(CHR = vcat(string.(1:22), ["X", "Y"]))
    if build == 38
        storage.maxpos = [GRCh38_totlength[chr] for chr in storage.CHR]
    else
        storage.maxpos = [GRCh37_totlength[chr] for chr in storage.CHR]
    end
    storage.add = cumsum(storage.maxpos) - storage.maxpos
    df = leftjoin(df, storage; on = :CHR)
    df.x = df.BP + df.add
    xmax = sum(unique(df.maxpos))
    indeven = findall(in(vcat(string.(1:2:22), "X")), df.CHR)
    indodd = findall(in(vcat(string.(2:2:23), "Y")), df.CHR)
    scatter!(ax, view(df, indeven, :x), view(df, indeven, :P), markersize = 1.5, color = chromcolors[1])
    scatter!(ax, view(df, indodd, :x), view(df, indodd, :P), markersize = 1.5, color = chromcolors[2])
    if !isnothing(scattercolor)
        ind = df.P .> -log(10, p)
        dfsig = view(df, ind, :)
        scatter!(ax, dfsig.x, dfsig.P, markersize = 1.5, color = scattercolor)
    end
    if !isnothing(linecolor)
        lines!(ax, [0.0, xmax], fill(-log(10, p), 2), color = linecolor, linewidth = 0.75)
    end
    xlims!(ax, 0, xmax)
    ylims!(ax, 0, ymax)
    hideydecorations!(ax, label = false, ticklabels = false, ticks = false)
    hidexdecorations!(ax, label = false, ticklabels = false)
    ax.xlabel = "Chromosome"
    ax.ylabel = rich("-log", subscript("10"), rich("P", font = :italic))
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