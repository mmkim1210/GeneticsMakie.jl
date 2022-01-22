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

function calcluateld!(
    gwas::DataFrame, 
    ref::SnpData; 
    snp::Union{AbstractString, Tuple{AbstractString, Int}} = "index"
)

    gwas.ind = findmissing(findsnps(gwas, ref))
    dropmissing!(gwas, "ind")
    n = size(gwas, 1)
    gwas.LD = fill(0.0, n)
    if snp == "index"
        i = argmax(gwas.P)
        snp = gwas.SNP[i]
        geno = convert(Matrix{Float64}, ref.snparray[:, gwas.ind])
        for j in 1:n
            gwas.LD[j] = cor(geno[:, i], geno[:, j])^2    
        end
        gwas.index = fill(snp, n)
        return
    else
        if snp isa AbstractString 
            ((chr, bp) = getsnpinfo(snp, ref))
        else 
            ((chr, bp) = snp)
        end
        i = findfirst((gwas.CHR .== chr) .& (gwas.BP .== bp))
        if isnothing(i)
            snp = ""
            i = findfirst((ref.snp_info.chromosome .== chr) .& (ref.snp_info.position .== bp))
            geno = convert(Matrix{Float64}, ref.snparray[:, [i; gwas.ind]])
            for j in 1:n
                gwas.LD[j] = cor(geno[:, 1], geno[:, j + 1])^2    
            end
        else
            snp = gwas.SNP[i]
            geno = convert(Matrix{Float64}, ref.snparray[:, gwas.ind])
            for j in 1:n
                gwas.LD[j] = cor(geno[:, i], geno[:, j])^2    
            end
        end
        gwas.index = fill(snp, n)
        return
    end
end

"""
    plotlocus!(ax::Axis, chromosome::AbstractString, range1::Real, range2::Real, gwas::DataFrame; kwargs)
    plotlocus!(ax::Axis, chromosome::AbstractString, bp::Real, gwas::DataFrame; kwargs)
    plotlocus!(ax::Axis, gene::AbstractString, gwas::DataFrame, gencode::DataFrame; kwargs)

Plot `gwas` results within a given `chromosome` and genomic range between `range1` 
and `range2`.

Alternatively, plot within a given `chromosome` and a certain `window` around a 
genomic coordinate `bp` or plot within a certain `window` around `gene`.

# Arguments
- `ld::Union{Nothing, SnpData, Tuple{SnpData, Union{AbstractString, Tuple{AbstractString, Int}}}} = nothing`: 
    the reference panel for which LD is calculated.
- `ymax::Real`: the maximum value for y axis. 
- `window::Real = 1e6`: the window around `bp` or `gene`. 
"""
function plotlocus!(
    ax::Axis,
    chromosome::AbstractString,
    range1::Real,
    range2::Real,
    gwas::DataFrame;
    ld::Union{Nothing, SnpData, Tuple{SnpData, Union{AbstractString, Tuple{AbstractString, Int}}}} = nothing,
    ymax::Real = 0
)

    df = filter(x -> (x.CHR == chromosome) && (x.BP >= range1) && (x.BP <= range2), gwas)
    if nrow(df) == 0
        ymax == 0 ? ymax = 10 : nothing
        yticks = setticks(ymax)
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
        return
    end
    df.P = -log.(10, df.P)
    if ymax == 0
        ymax = maximum(df.P) / 4 * 5
        ymax <= 10 ? ymax = 10 : nothing 
        yticks = setticks(ymax)
    else
        yticks = setticks(ymax)
    end
    if !isnothing(ld)
        typeof(ld) == SnpData ? calcluateld!(df, ld) : calcluateld!(df, ld[1]; snp = ld[2])
        scatter!(ax, df.BP, df.P, color = df.LD, colorrange = (0, 1),
            colormap = (:gray60, :red2), markersize = 1.5)
        if typeof(ld) == SnpData
            ind = argmax(df.P)
            bp = df.BP[ind]
            p = df.P[ind]
            scatter!(ax, [bp], [p], color = :purple1, markersize = 4.0, marker = '◆')
            text!(ax, "$(df.index[1])", position = (bp, p), textsize = 6, align = (:center, :bottom))    
        elseif length(df.index[1]) > 0
            ind = findfirst(df.SNP .== df.index[1])
            bp, p = df.BP[ind], df.P[ind]
            scatter!(ax, [bp], [p], color = :purple1, markersize = 4.0, marker = '◆')
            text!(ax, "$(df.index[1])", position = (bp, p), textsize = 6, align = (:center, :bottom))    
        end
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

plotlocus!(ax::Axis, chromosome::AbstractString, bp::Real, gwas::DataFrame; window::Real = 1e6, kwargs...) =
    plotlocus!(ax, chromosome, bp - window, bp + window, gwas; kwargs...)

function plotlocus!(ax::Axis, gene::AbstractString, gwas::DataFrame, gencode::DataFrame; window::Real = 1e6, kwargs...)
    chr, start, stop = findgene(gene, gencode)
    plotlocus!(ax, chr, start - window, stop + window, gwas; kwargs...)
end