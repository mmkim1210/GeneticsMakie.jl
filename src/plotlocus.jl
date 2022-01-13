function setticks(y::Real)
    if y == 10
        return 0:3:10
    elseif 10 < y <= 20
        return 0:5:y
    elseif 20 < y <= 50
        return 0:10:y
    else
        s = div(y, 4)
        d, r = divrem(s, 10)
        return r < 5 ? (0:10d:y) : (0:((d + 1) * 10):y)
    end
end

"""
    calcluateld!(gwas::DataFrame, ref::SnpData; snp)

Calculate LD between the most significant SNP and other SNPs in `gwas` by using
`ref`. Optionally, the target SNP can be switched to `snp`.
"""
function calcluateld!(gwas::DataFrame, ref::SnpData; snp::AbstractString = "index")
    gwas.ind = findmissing(findsnps(gwas, ref))
    dropmissing!(gwas, "ind")
    n = size(gwas, 1)
    gwas.LD = fill(0.0, n)
    if snp == "index"
        i = argmax(gwas.P)
        snp = gwas.SNP[i]
    else
        i = findfirst(gwas.SNP .== snp)
    end
    gwas.index = fill(snp, n)
    geno = convert(Matrix{Float64}, ref.snparray[:, gwas.ind])
    for j in 1:n
        gwas.LD[j] = cor(geno[:, i], geno[:, j])^2    
    end
    return
end

"""
    plotlocus!(ax::Axis, chromosome::AbstractString, range1::Real, range2::Real, gwas::DataFrame; colorld, ref, snp, ymax)

Plot `gwas` results within a given `chromosome` and genomic range between `range1` 
and `range2`. Optionally, SNPs can be colored by LD via `colorld` using `ref`.
The default SNP for which LD is calculated is index SNP, which can be changed to `snp`.
"""
function plotlocus!(ax::Axis,
    chromosome::AbstractString,
    range1::Real,
    range2::Real,
    gwas::DataFrame;
    colorld::Bool = false,
    ref::Union{Nothing, SnpData} = nothing,
    snp::Union{AbstractString, Tuple{Int, Int}} = "index",
    ymax::Real = 0)

    df = filter(x -> (x.CHR == chromosome) && (x.BP >= range1) && (x.BP <= range2), gwas)
    df.P = -log.(10, df.P)
    if ymax == 0
        ymax = maximum(df.P) / 4 * 5
        ymax <= 10 ? ymax = 10 : nothing 
        yticks = setticks(ymax)
    else
        yticks = setticks(ymax)
    end
    if colorld
        snp == "index" ? calcluateld!(df, ref) : calcluateld!(df, ref; snp = snp)
        scatter!(ax, df.BP, df.P, color = df.LD, colorrange = (0, 1),
            colormap = (:gray60, :red2), markersize = 1.5)
        if snp == "index"
            bp = getindex(df.BP, argmax(df.P))
            p = maximum(df.P)
        else
            ind = findfirst(df.SNP .== snp)
            bp, p = df.BP[ind], df.P[ind]    
        end
        scatter!(ax, [bp], [p], color = :purple1, markersize = 4.0, marker = '◆')
        text!(ax, "$(df.index[1])", position = (bp, p), textsize = 6, align = (:center, :bottom))
    else
        scatter!(ax, df.BP, df.P, color = :gray60, markersize = 1.5)
    end
    ax.spinewidth = 0.75
    ax.ytickwidth = 0.75
    ax.ylabelsize = 6
    ax.yticklabelsize = 6
    ax.yticksize = 3
    ax.yticks = setticks(ymax)
    xlims!(ax, range1, range2)
    ylims!(ax, 0, ymax)
    hidespines!(ax, :t, :r)
    hidexdecorations!(ax)
    hideydecorations!(ax, ticks = false, label = false, ticklabels = false)
end

"""
    plotlocus!(ax::Axis, chromosome::AbstractString, bp::Real, gwas::DataFrame; window, colorld, ref, snp, ymax, title)

Plot `gwas` results within a given `chromosome` and a certain `window` around a 
genomic coordinate `bp`. The default window is 1 Mb.
"""
plotlocus!(ax::Axis, chromosome::AbstractString, bp::Real, gwas::DataFrame; window::Real = 1e6, kwargs...) =
    plotlocus!(ax, chromosome, bp - window, bp + window, gwas; kwargs...)

"""
    plotlocus!(ax::Axis, gene::AbstractString, gwas::DataFrame, gencode::DataFrame; window, colorld, ref, snp, ymax, title)

Plot `gwas` results within a certain window around `gene`. The default window is 1 Mb.
"""
function plotlocus!(ax::Axis, gene::AbstractString, gwas::DataFrame, gencode::DataFrame; window::Real = 1e6, kwargs...)
    ind = findfirst(isequal(gene), gencode.gene_name)
    plotlocus!(ax, gencode.seqnames[ind], gencode.start[ind] - window, gencode[ind, :end] + window, gwas; kwargs...)
end    