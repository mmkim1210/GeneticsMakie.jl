findlocus(CHR::AbstractVector, BP::AbstractVector, chr::AbstractString, range1::Real, range2::Real) =
    ind = findall((CHR .== chr) .& (BP .>= range1) .& (BP .<= range2))

findlocus(gwas::DataFrame, chr::AbstractString, range1::Real, range2::Real) =
    findlocus(gwas.CHR, gwas.BP, chr, range1, range2)

findlocus(ref::SnpData, chr::AbstractString, range1::Real, range2::Real) =
    findlocus(ref.snp_info.chromosome, ref.snp_info.position, chr, range1, range2)

function mungenames!(gwas::DataFrame)
    rename!(uppercase, gwas)
    colnames = Dict(
        "SNPID" => "SNP",
        "MARKER" => "SNP",
        "MARKERNAME" => "SNP",
        "RSID" => "SNP",
        "ID" => "SNP",
        "CHROMOSOME" => "CHR",
        "CHROM" => "CHR",
        "POSITION" => "BP",
        "POS" => "BP",
        "ALLELE1" => "A1",
        "EFFECT_ALLELE" => "A1",
        "ALLELE1" => "A1",
        "ALLELE2" => "A2",
        "OTHER_ALLELE" => "A2",
        "NON_EFFECT_ALLELE" => "A2",
        "ZSCORE" => "Z",
        "PVAL" => "P",
        "PVALUE" => "P"
    )
    for name in names(gwas)
        haskey(colnames, name) ? rename!(gwas, name => colnames[name]) : nothing
    end
end

function mungetypes!(gwas::DataFrame)
    eltype(gwas.CHR) <: AbstractString ? nothing : gwas.CHR = string.(gwas.CHR)
end

function mungesumstats!(gwas::Vector{DataFrame})
    for i in eachindex(gwas)
        mungenames!(gwas[i])
        mungetypes!(gwas[i])
    end
end

mungesumstats!(gwas::DataFrame) = mungesumstats!([gwas])

function matchsnps(CHR::AbstractVector, BP::AbstractVector, chr::AbstractString, bp::Real, window)
    ind = Int[]
    for i in eachindex(CHR)
        if CHR[i] == chr && BP[i] < bp + window && BP[i] > bp - window
            push!(ind, i)
        end
    end
    return ind
end

function findgwasloci(gwas::DataFrame; p = 5e-8)
    loci = DataFrame(CHR = String[], BP = Int[], P = Float64[])
    df = gwas[findall(x -> x < p, gwas.P), ["CHR", "BP", "P"]]
    while nrow(df) > 0
        ind = argmin(df.P)
        push!(loci, df[ind, :])
        ind = matchsnps(df.CHR, df.BP, df.CHR[ind], df.BP[ind], 1e6)
        isnothing(ind) ? nothing : deleteat!(df, ind)
    end
    return loci
end

function findgwasloci(gwas::Vector{DataFrame}; kwargs...)
    loci = Vector{DataFrame}(undef, length(gwas))
    for i in eachindex(gwas)
        loci[i] = findgwasloci(gwas[i]; kwargs...)
    end
    for i in 2:length(gwas)
        for j in 1:nrow(loci[1])
            ind = matchsnps(loci[i].CHR, loci[i].BP, loci[1].CHR[j], loci[1].BP[j], 1e6)
            isnothing(ind) ? nothing : deleteat!(loci[i], ind)
        end
    end
    return vcat(loci...)
end

function findsnps(gwas::DataFrame, ref::SnpData)
    ind = Vector{Union{Missing, Int}}(undef, nrow(gwas))
    for i in eachindex(ind)
        j = findfirst((ref.snp_info.chromosome .== gwas.CHR[i]) .& (ref.snp_info.position .== gwas.BP[i]))
        isnothing(j) ? ind[i] == missing : ind[i] = j
    end
    return ind
end