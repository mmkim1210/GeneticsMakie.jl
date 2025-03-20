function mungenames!(gwas::DataFrame)
    rename!(uppercase, gwas)
    colnames = Dict(
        "SNPID" => "SNP",
        "MARKER" => "SNP",
        "MARKERNAME" => "SNP",
        "RSID" => "SNP",
        "RS_ID" => "SNP",
        "RS_NUMBER" => "SNP",
        "RS_NUMBERS" => "SNP",
        "ID" => "SNP",
        "VARIANT_ID" => "SNP",
        "#CHR" => "CHR",
        "#CHROM" => "CHR",
        "CHROMOSOME" => "CHR",
        "CHROM" => "CHR",
        "HG19CHR" => "CHR",
        "CHR.ICOGS" => "CHR",
        "POSITION" => "BP",
        "POS" => "BP",
        "POSITION(HG19)" => "BP",
        "BP_HG19" => "BP",
        "BASE_PAIR_LOCATION" => "BP",
        "POS_B37" => "BP",
        "POSGRCH37" => "BP",
        "POSITION.ICOGS" => "BP",
        "ALLELE1" => "A1",
        "ALLELE_1" => "A1",
        "EFFECT_ALLELE" => "A1",
        "TESTED_ALLELE" => "A1",
        "TESTEDALLELE" => "A1",
        "ALLELE1" => "A1",
        "REF" => "A1",
        "REFERENCE_ALLELE" => "A1",
        "EA" => "A1",
        "INC_ALLELE" => "A1",
        "EFFECT.META" => "A1",
        "ALLELE2" => "A2",
        "ALLELE_2" => "A2",
        "OTHER_ALLELE" => "A2",
        "OTHERALLELE" => "A2",
        "NON_EFFECT_ALLELE" => "A2",
        "NONEFFECT_ALLELE" => "A2",
        "ALT" => "A2",
        "NEA" => "A2",
        "BASELINE.META" => "A2",
        "DEC_ALLELE" => "A2",
        "ZSCORE" => "Z",
        "GC_ZSCORE" => "Z",
        "Z-SCORE" => "Z",
        "T-STATISTIC" => "Z",
        "PVAL" => "P",
        "PVALUE" => "P",
        "P-VALUE" => "P",
        "P-VAL" => "P",
        "P_DGC" => "P",
        "P_VALUE" => "P",
        "GC.PVALUE" => "P",
        "GC_PVALUE" => "P",
        "P.VALUE" => "P",
        "P.2GC" => "P",
        "ALL_INV_VAR_META_P" => "P",
        "P.META" => "P",
        "OR(A1)" => "OR",
        "OR(MINALLELE)" => "OR",
        "ORX" => "OR",
        "SE_DGC" => "SE",
        "STDERR" => "SE",
        "STANDARD_ERROR" => "SE",
        "ALL_INV_VAR_META_SEBETA" => "SE",
        "STDERRLOGOR" => "SE",
        "SE.2GC" => "SE",
        "SDE.META" => "SE",
        "EFFECT" => "BETA",
        "EFFECTS" => "BETA",
        "STDBETA" => "BETA",
        "B" => "BETA",
        "LOGOR" => "BETA",
        "LOG_ODDS" => "BETA",
        "ALL_INV_VAR_META_BETA" => "BETA",
        "BETA.META" => "BETA",
        "EAF" => "FRQ",
        "AF" => "FRQ",
        "FREQ1" => "FRQ",
        "FREQ" => "FRQ",
        "EAF_A1" => "FRQ",
        "EAF_UKB" => "FRQ",
        "MAF" => "FRQ",
        "EFFECT_ALLELE_FREQ" => "FRQ",
        "EFFECT_ALLELE_FREQUENCY" => "FRQ",
        "FREQ_TESTED_ALLELE_IN_HRS" => "FRQ",
        "EAF_EUR_UKB" => "FRQ",
        "FREQ_HAPMAP" => "FRQ",
        "EAF_HRC" => "FRQ",
        "NUMBER_OF_STUDIES" => "NSTUDY",
        "NSTUDIES" => "NSTUDY",
        "N_STUDIES" => "NSTUDY",
        "N_STUDY" => "NSTUDY",
        "NCASES" => "N_CAS",
        "NCASE" => "N_CAS",
        "N_CASE" => "N_CAS",
        "NCAS" => "N_CAS",
        "NCA" => "N_CAS",
        "N_CASES" => "N_CAS",
        "CASES_N" => "N_CAS",
        "NCONTROLS" => "N_CON",
        "NCONTROL" => "N_CON",
        "N_CONTROLS" => "N_CON",
        "CONTROLS_N" => "N_CON",
        "NCON" => "N_CON",
        "NCO" => "N_CON",
        "N_CONTROL" => "N_CON",
        "ALL_META_N" => "N",
        "N_ANALYZED" => "N",
        "IMPINFO" => "INFO",
        "INFO_UKB" => "INFO",
        "MEDIAN_INFO" => "INFO",
        "MININFO" => "INFO"
    )
    if in("SNP", [get(colnames, name, 0) for name in names(gwas)]) && "SNP" in names(gwas)
        rename!(gwas, "SNP" => "EXTRAID")
    end
    for name in names(gwas)
        haskey(colnames, name) ? rename!(gwas, name => colnames[name]) : nothing
    end
end

function mungesnpid!(gwas::DataFrame)
    if !("SNP" in names(gwas))
        gwas.SNP = string.(gwas.CHR, ":", gwas.BP)
    end
    if any(ismissing.(gwas.SNP))
        ind = ismissing.(gwas.SNP)
        gwas.SNP = string.(gwas.SNP)
        gwas.SNP[ind] = string.(gwas.CHR[ind], ":", gwas.BP[ind])
    end
end

function mungetypes!(gwas::DataFrame)
    eltype(gwas.CHR) <: AbstractString ? nothing : gwas.CHR = string.(gwas.CHR)
    if eltype(gwas.BP) != Int && eltype(gwas.BP) <: AbstractString
        gwas.BP = parse.(Int, gwas.BP)
    elseif eltype(gwas.BP) != Int && eltype(gwas.BP) <: Real
        gwas.BP = convert.(Int, gwas.BP)
    end
end

function mungealleles!(gwas::DataFrame)
    gwas.A1 = uppercase.(gwas.A1)
    gwas.A2 = uppercase.(gwas.A2)
end

function mungezscore!(gwas::DataFrame)
    if "BETA" in names(gwas) && "SE" in names(gwas)
        gwas.Z = gwas.BETA ./ gwas.SE
        filter!(x -> !isnan(x.Z), gwas)
    elseif "OR" in names(gwas) && "SE" in names(gwas)
        gwas.Z = log.(gwas.OR) ./ gwas.SE
        filter!(x -> !isnan(x.Z), gwas)
    elseif "BETA" in names(gwas)
        gwas.Z .= 0.0
        for i in 1:nrow(gwas)
            if gwas.BETA[i] >= 0
                gwas.Z[i] = sqrt(cquantile(Chisq(1), gwas.P[i]))
            else
                gwas.Z[i] = -sqrt(cquantile(Chisq(1), gwas.P[i]))
            end
        end
    elseif "OR" in names(gwas)
        gwas.BETA = log.(gwas.OR)
        gwas.Z .= 0.0
        for i in 1:nrow(gwas)
            if gwas.BETA[i] >= 0
                gwas.Z[i] = sqrt(cquantile(Chisq(1), gwas.P[i]))
            else
                gwas.Z[i] = -sqrt(cquantile(Chisq(1), gwas.P[i]))
            end
        end
    elseif "OVERALL" in names(gwas)
        gwas.Z .= 0.0
        for i in 1:nrow(gwas)
            if gwas.OVERALL[i] == "+"
                gwas.Z[i] = sqrt(cquantile(Chisq(1), gwas.P[i]))
            else
                gwas.Z[i] = -sqrt(cquantile(Chisq(1), gwas.P[i]))
            end
        end
    end
end

function mungepvalue!(gwas::DataFrame)
    filter!(x -> 0 <= x.P <= 1, gwas)
    gwas.P = clamp.(gwas.P, floatmin(0.0), 1)
end

"""
    mungesumstats!(gwas::DataFrame)
    mungesumstats!(gwas::Vector{DataFrame})

Munge `gwas` by harmonizing the names of columns, their types, and P values, among others.
"""
function mungesumstats!(gwas::Vector{DataFrame})
    for i in eachindex(gwas)
        mungenames!(gwas[i])
        ind = findall(in(["SNP", "CHR", "BP", "A1", "A2", "Z", "P", "BETA", "OR", "SE", "OVERALL"]), names(gwas[i]))
        select!(gwas[i], ind)
        mungesnpid!(gwas[i])
        dropmissing!(gwas[i])
        mungetypes!(gwas[i])
        if "A1" in names(gwas[i]) && "A2" in names(gwas[i])
            mungealleles!(gwas[i])
        end
        mungezscore!(gwas[i])
        mungepvalue!(gwas[i])
        if "A1" in names(gwas[i]) && "A2" in names(gwas[i]) && "Z" in names(gwas[i])
            select!(gwas[i], :SNP, :CHR, :BP, :A1, :A2, :Z, :P)
        elseif "A1" in names(gwas[i]) && "A2" in names(gwas[i])
            select!(gwas[i], :SNP, :CHR, :BP, :A1, :A2, :P)
        elseif "Z" in names(gwas[i])
            select!(gwas[i], :SNP, :CHR, :BP, :Z, :P)
        else
            select!(gwas[i], :SNP, :CHR, :BP, :P)
        end
    end
end

mungesumstats!(gwas::DataFrame) = mungesumstats!([gwas])

findlocus(CHR::AbstractVector, BP::AbstractVector, chr::AbstractString, range1::Real, range2::Real) =
    findall((CHR .== chr) .& (BP .>= range1) .& (BP .<= range2))

findlocus(gwas::DataFrame, chr::AbstractString, range1::Real, range2::Real) =
    findlocus(gwas.CHR, gwas.BP, chr, range1, range2)

findlocus(ref::SnpData, chr::AbstractString, range1::Real, range2::Real) =
    findlocus(ref.snp_info.chromosome, ref.snp_info.position, chr, range1, range2)

"""
    findgwasloci(gwas::DataFrame)
    findgwasloci(gwas::Vector{DataFrame})

Find genome-wide significant loci for `gwas` that are separated from each
other by at least 1 Mb.

Alternatively, find genome-wide significant loci across multiple `gwas` that
are all separated by at least 1 Mb.

# Keyword argument
```
p::Real        genome-wide significance threshold; default 5e-8
```
"""
function findgwasloci(gwas::DataFrame; p::Real = 5e-8)
    loci = DataFrame(CHR = String[], BP = Int[], P = Float64[])
    df = gwas[findall(x -> x < p, gwas.P), [:CHR, :BP, :P]]
    while nrow(df) > 0
        ind = argmin(df.P)
        push!(loci, df[ind, :])
        ind = findlocus(df.CHR, df.BP, df.CHR[ind], df.BP[ind] - 1e6, df.BP[ind] + 1e6)
        isnothing(ind) ? nothing : deleteat!(df, ind)
    end
    return loci
end

function findgwasloci(gwas::Vector{DataFrame}; kwargs...)
    loci = Vector{DataFrame}(undef, length(gwas))
    for i in eachindex(gwas)
        loci[i] = findgwasloci(gwas[i]; kwargs...)
    end
    for i in 1:length(gwas)
        for j in (i + 1):length(gwas)
            for k in 1:nrow(loci[i])
                ind = findlocus(loci[j].CHR, loci[j].BP, loci[i].CHR[k], loci[i].BP[k] - 1e6, loci[i].BP[k] + 1e6)
                isnothing(ind) ? nothing : deleteat!(loci[j], ind)
            end
        end
    end
    return vcat(loci...)
end

function findsnps(CHR₁::AbstractVector, BP₁::AbstractVector, CHR₂::AbstractVector, BP₂::AbstractVector)
    ind = Vector{Union{Missing, Int}}(undef, length(CHR₁))
    for i in 1:length(ind)
        j = findfirst((CHR₂ .== CHR₁[i]) .& (BP₂ .== BP₁[i]))
        isnothing(j) ? ind[i] = missing : ind[i] = j
    end
    return ind
end

function findsnps(
    CHR₁::AbstractVector,
    BP₁::AbstractVector,
    A1₁::AbstractVector,
    A2₁::AbstractVector,
    CHR₂::AbstractVector,
    BP₂::AbstractVector,
    A1₂::AbstractVector,
    A2₂::AbstractVector;
    matchalleles::Bool = true
    )
    if matchalleles
        ind = Matrix{Union{Missing, Int}}(undef, length(CHR₁), 2)
        for i in 1:size(ind, 1)
            j = findfirst((CHR₂ .== CHR₁[i]) .& (BP₂ .== BP₁[i]))
            if isnothing(j)
                ind[i, :] .= missing
            else
                if A1₂[j] == A1₁[i] && A2₂[j] == A2₁[i]
                    ind[i, 1], ind[i, 2] = j, missing
                elseif A2₂[j] == A1₁[i] && A1₂[j] == A2₁[i]
                    ind[i, 1], ind[i, 2] = missing, j
                else
                    ind[i, :] .= missing
                end
            end
        end
        return ind
    else
        return findsnps(CHR₁, BP₁, CHR₂, BP₂)
    end
end

function findsnps(gwas::DataFrame, ref::SnpData; matchalleles::Bool = true)
    if matchalleles
        return findsnps(gwas.CHR, gwas.BP, gwas.A1, gwas.A2,
            ref.snp_info.chromosome, ref.snp_info.position, ref.snp_info.allele1, ref.snp_info.allele2)
    else
        return findsnps(gwas.CHR, gwas.BP, ref.snp_info.chromosome, ref.snp_info.position)
    end
end

function findsnps(ref::SnpData, gwas::DataFrame; matchalleles::Bool = true)
    if matchalleles
        return findsnps(ref.snp_info.chromosome, ref.snp_info.position, ref.snp_info.allele1, ref.snp_info.allele2,
            gwas.CHR, gwas.BP, gwas.A1, gwas.A2)
    else
        return findsnps(ref.snp_info.chromosome, ref.snp_info.position, gwas.CHR, gwas.BP)
    end
end

function findsnps(gwas₁::DataFrame, gwas₂::DataFrame; matchalleles::Bool = true)
    if matchalleles
        return findsnps(gwas₁.CHR, gwas₁.BP, gwas₁.A1, gwas₁.A2, gwas₂.CHR, gwas₂.BP, gwas₂.A1, gwas₂.A2)
    else
        return findsnps(gwas₁.CHR, gwas₁.BP, gwas₂.CHR, gwas₂.BP)
    end
end

"""
    findclosestgene(chr::AbstractString, bp::Real, gencode::DataFrame)
    findclosestgene(df::DataFrame, gencode::DataFrame)

Find the closest gene to a genomic coordinate or a list of genomic coordinates using `gencode`.

# Keyword arguments
```
start::Bool             find closest gene from gene start site; default false
proteincoding::Bool     find closest protein coding gene; default false
```
"""
function findclosestgene(
    chr::AbstractString,
    bp::Real,
    gencode::DataFrame;
    start::Bool = false,
    proteincoding::Bool = false
    )
    if proteincoding
        df = filter(x -> (x.seqnames == chr) && (x.feature == "gene") && (x.gene_type == "protein_coding"), gencode)
    else
        df = filter(x -> (x.seqnames == chr) && (x.feature == "gene"), gencode)
    end
    if start
        df.dist₁ .= 1
        for i in 1:nrow(df)
            df.strand[i] == "+" ? df.dist₁[i] = abs(df.start[i] - bp) : df.dist₁[i] = abs(df.end[i] - bp)
        end
        ind₁ = argmin(df.dist₁)
        return df.gene_name[ind₁], df.dist₁[ind₁]
    else
        df.dist₁ = abs.(df.start .- bp)
        df.dist₂ = abs.(df.end .- bp)
        ind₁ = argmin(df.dist₁)
        ind₂ = argmin(df.dist₂)
        if df.dist₁[ind₁] < df.dist₂[ind₂]
            return df.gene_name[ind₁], df.dist₁[ind₁]
        elseif df.dist₁[ind₁] > df.dist₂[ind₂]
            return df.gene_name[ind₂], df.dist₂[ind₂]
        elseif ind₁ == ind₂
            return df.gene_name[ind₁], df.dist₁[ind₁]
        elseif df.strand[ind₁] == "+" && df.strand[ind₂] == "+"
            return df.gene_name[ind₁], df.dist₁[ind₁]
        elseif df.strand[ind₁] == "-" && df.strand[ind₂] == "-"
            return df.gene_name[ind₂], df.dist₁[ind₂]
        elseif df.gene_type[ind₁] == "protein_coding" && df.gene_type[ind₂] != "protein_coding"
            return df.gene_name[ind₁], df.dist₁[ind₁]
        elseif df.gene_type[ind₁] != "protein_coding" && df.gene_type[ind₂] == "protein_coding"
            return df.gene_name[ind₂], df.dist₁[ind₂]
        else
            return df.gene_name[ind₁], df.dist₁[ind₁]
        end
    end
end

function findclosestgene(CHR::AbstractVector, BP::AbstractVector, gencode::DataFrame; kwargs...)
    gencodeₛ = select(gencode, :seqnames, :start, :end, :strand, :gene_type, :feature, :gene_name)
    filter!(x -> x.feature == "gene", gencodeₛ)
    storage = DataFrame(CHR = String[], BP = Int[], gene = String[], distance = Int[])
    for i in eachindex(CHR)
        gene, dist = findclosestgene(CHR[i], BP[i], gencodeₛ; kwargs...)
        push!(storage, [CHR[i], BP[i], gene, dist])
    end
    return storage
end

findclosestgene(df::DataFrame, gencode::DataFrame; kwargs...) = findclosestgene(df.CHR, df.BP, gencode; kwargs...)

function findclosestgenes(
    chr::AbstractString,
    bp::Real,
    gencode::DataFrame;
    start::Bool = false,
    proteincoding::Bool = false,
    n::Real = 5
    )
    if proteincoding
        df = filter(x -> (x.seqnames == chr) && (x.feature == "gene") && (x.gene_type == "protein_coding"), gencode)
    else
        df = filter(x -> (x.seqnames == chr) && (x.feature == "gene"), gencode)
    end
    if start
        df.dist₁ .= 1
        for i in 1:nrow(df)
            df.strand[i] == "+" ? df.dist₁[i] = abs(df.start[i] - bp) : df.dist₁[i] = abs(df.end[i] - bp)
        end
        sort!(df, :dist₁)
        select!(df, :gene_name, :dist₁)
        rename!(df, :dist₁ => :distance)
        return first(df, n)
    else
        df.dist₁ = abs.(df.start .- bp)
        df.dist₂ = abs.(df.end .- bp)
        df.dist = min.(df.dist₁, df.dist₂)
        sort!(df, :dist)
        select!(df, :gene_name, :dist)
        rename!(df, :dist => :distance)
        return first(df, n)
    end
end

function getsnpinfo(snp::AbstractString, SNP::AbstractVector, CHR::AbstractVector, BP::AbstractVector)
    ind = findfirst(isequal(snp), SNP)
    isnothing(ind) ? nothing : (CHR[ind], BP[ind])
end

getsnpinfo(snp::AbstractString, df::DataFrame) = getsnpinfo(snp, df.SNP, df.CHR, df.BP)

getsnpinfo(snp::AbstractString, ref::SnpData) =
    getsnpinfo(snp, ref.snp_info.snpid, ref.snp_info.chromosome, ref.snp_info.position)

function getsnpinfo(chr::AbstractString, bp::Real, SNP::AbstractVector, CHR::AbstractVector, BP::AbstractVector)
    ind = findfirst((CHR .== chr) .& (BP .== bp))
    isnothing(ind) ? nothing : SNP[ind]
end

getsnpinfo(chr::AbstractString, bp::Real, df::DataFrame) = getsnpinfo(chr, bp, df.SNP, df.CHR, df.BP)

getsnpinfo(chr::AbstractString, bp::Real, ref::SnpData) =
    getsnpinfo(chr, bp, ref.snp_info.snpid, ref.snp_info.chromosome, ref.snp_info.position)

function flipalleles(gwas::Vector{DataFrame}; outer = false)
    if outer
        df = outerjoin(gwas..., on = [:CHR, :BP], makeunique = true)
        for i in 1:(length(gwas) - 1)
            for j in 1:nrow(df)
                if ismissing(df.A1[j]) || ismissing(df.A2[j]) || ismissing(df[j, "A1_$i"]) || ismissing(df[j, "A2_$i"])
                    continue
                elseif df.A1[j] == df[j, "A1_$i"] && df.A2[j] == df[j, "A2_$i"]
                    continue
                elseif df.A1[j] == df[j, "A2_$i"] && df.A2[j] == df[j, "A1_$i"]
                    df[j, "Z_$i"] = -df[j, "Z_$i"]
                else
                    df[j, "Z_$i"] = NaN
                end
            end
        end
    else
        df = innerjoin(gwas..., on = [:CHR, :BP], makeunique = true)
        for i in 1:(length(gwas) - 1)
            for j in 1:nrow(df)
                if df.A1[j] == df[j, "A1_$i"] && df.A2[j] == df[j, "A2_$i"]
                    continue
                elseif df.A1[j] == df[j, "A2_$i"] && df.A2[j] == df[j, "A1_$i"]
                    df[j, "Z_$i"] = -df[j, "Z_$i"]
                else
                    df[j, "Z_$i"] = NaN
                end
            end
        end
    end
    return df
end

flipalleles(gwas₁::DataFrame, gwas₂::DataFrame; kwargs...) = flipalleles([gwas₁, gwas₂]; kwargs...)
