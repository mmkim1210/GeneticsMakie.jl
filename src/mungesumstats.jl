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
        gwas.Z = fill(0.0, nrow(gwas))
        for i in 1:nrow(gwas)
            if gwas.BETA[i] >= 0
                gwas.Z[i] = sqrt(cquantile(Chisq(1), gwas.P[i]))
            else
                gwas.Z[i] = -sqrt(cquantile(Chisq(1), gwas.P[i]))
            end
        end
    elseif "OR" in names(gwas)
        gwas.BETA = log.(gwas.OR)
        gwas.Z = fill(0.0, nrow(gwas))
        for i in 1:nrow(gwas)
            if gwas.BETA[i] >= 0
                gwas.Z[i] = sqrt(cquantile(Chisq(1), gwas.P[i]))
            else
                gwas.Z[i] = -sqrt(cquantile(Chisq(1), gwas.P[i]))
            end
        end
    elseif "OVERALL" in names(gwas)
        gwas.Z = fill(0.0, nrow(gwas))
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
    findgwasloci(gwas::DataFrame; p::Real)
    findgwasloci(gwas::Vector{DataFrame}; p::Real)

Find genome-wide significant loci for `gwas` that are separated from each
other by at least 1 Mb.

Alternatively, find genome-wide significant loci across multiple `gwas` that 
are all separated by at least 1 Mb. `p` determines the genome-wide significance threshold, 
which is 5e-8 by default.
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
    findclosestgene(chr::AbstractString, bp::Real, gencode::DataFrame; start::Bool, proteincoding::Bool)
    findclosestgene(df::DataFrame, gencode::DataFrame; start::Bool, proteincoding::Bool)

Find the closest gene(s) to a genomic coordinate or a list of genomic coordinates using `gencode`. 

Optionally, the closest gene can be defined from the gene start site using `start`,
and only protein coding genes can be considered using `proteincoding`. 
The default `start` and `proteincoding` are `false`.
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

# Liftover summary statistics

function parsechain(chainpath)
    chains = NamedTuple{(:score, :tseq, :tsize, :tstrand, :tstart, :tend, :qseq, :qsize, 
                         :qstrand, :qstart, :qend, :alignments),
                        Tuple{Int64, String, Int64, String, Int64, Int64, String, 
                              Int64, String, Int64, Int64, Vector{Vector{Int64}}}}[]
    header = nothing
    alignments = Vector{Int64}[]
    magicnumbers = UInt8[]
    open(chainpath) do io
        readbytes!(io, magicnumbers, 2)
    end
    open(chainpath) do io
        # GZip decompression if the file is gzipped
        for line in eachline(magicnumbers == [0x1f, 0x8b] ? GzipDecompressorStream(io) : io)
            if line == ""
                continue
            end
            if startswith(line, "chain")
                linearr = split(line, r"[ \t]+")
                header = [parse(Int64, linearr[2]), # score
                          replace(string(linearr[3]), r"^chr" => ""), # tseq
                          parse(Int64, linearr[4]), # tsize
                          string(linearr[5]), # tstrand
                          parse(Int64, linearr[6]), # tstart
                          parse(Int64, linearr[7]), # tend
                          replace(string(linearr[8]), r"^chr" => ""), #qseq
                          parse(Int64, linearr[9]), # qsize
                          string(linearr[10]), # qstrand
                          parse(Int64, linearr[11]), # qstart
                          parse(Int64, linearr[12])] # qend
                continue
            end
            linearr = split(line, r"[ \t]+")
            push!(alignments, parse.(Int64, linearr))
            if length(linearr) == 1
                push!(chains,
                      (; zip([:score, :tseq, :tsize, :tstrand, :tstart, :tend, :qseq, :qsize, 
                              :qstrand, :qstart, :qend, :alignments],
                             [header..., alignments])...))
                alignments = Vector{Int64}[]
            end
        end
    end
    DataFrame(chains)
end

function expandchaindf(chaindf)
    vcat([let
              tpos = row.tstart + 1
              qpos = row.qstart + 1
              [let
                   tend = tpos + alignment[1] - 1
                   qend = qpos + alignment[1] - 1
                   newrow = (score = row.score,
                             tseq = row.tseq,
                             tsize = row.tsize,
                             tstrand = row.tstrand,
                             tstart = tpos,
                             tend = tend,
                             qseq = row.qseq,
                             qsize = row.qsize,
                             qstrand = row.qstrand,
                             qstart = qpos,
                             qend = qend)
                   if length(alignment) == 1
                       @assert tend == row.tend
                       @assert qend == row.qend
                   end
                   if length(alignment) == 3
                       tpos = tpos + alignment[1] + alignment[2]
                       qpos = qpos + alignment[1] + alignment[3]
                   end
                   newrow
               end
               for alignment in row.alignments]
          end
          for row in eachrow(chaindf)]...) |>
    DataFrame
end

"""
    readchain(path::AbstractString)

Read a chain file describing the genomic positions mapped between two 
reference genomes. Returns a DataFrame necessary for liftover. Sequence 
names are of type `String` and are stripped of the prefix "chr".
Information about the chain file format is present at 
[UCSC Genome Browser: Chain Format](https://genome.ucsc.edu/goldenPath/help/chain.html).
"""
function readchain(path::AbstractString)
    parsechain(path) |> expandchaindf |> df -> sort!(df, [:tseq, :tstart, :tend])
end

function findnewcoord(
        CHR::Vector{<:AbstractString}, BP::Vector{Int64}, chain::AbstractDataFrame;
        # Behavior when multiple positions map after liftover
        multiplematches::Symbol = :error # :error, :warning, :silent
    )
    @assert multiplematches in [:error, :warning, :silent]
    # Sort the positions so that we can iterate and liftover more efficiently
    chrbp = DataFrame(CHR = CHR, BP = BP)
    sortedidx = sortperm(chrbp, [:CHR, :BP])
    newpos = [NamedTuple{(:CHR, :BP, :strand, :score), Tuple{String, Int64, String, Int64}}[] for _ in 1:nrow(chrbp)]
    chainidxold = 1
    chainidxcur = 1
    chainidxend = nrow(chain)
    for idx in sortedidx
        chainidxcur = chainidxold
        chainidxset = false
        CHR::String, BP::Int64 = chrbp[idx, [:CHR, :BP]]
        while (chainidxcur <= chainidxend) && (chain[chainidxcur, :tseq] <= CHR)
            if CHR == chain[chainidxcur, :tseq]
                if chain[chainidxcur, :tstart] <= BP <= chain[chainidxcur, :tend]
                    # Update where to start looking in the chain file only once per 
                    # position
                    # A position might have multiple matches, and we don't want to 
                    # miss checking the regions that matched first
                    if !chainidxset
                        chainidxold = chainidxcur
                        chainidxset = true
                    end
                    # TODO: Handle "-" strands as well
                    @assert chain[chainidxcur, :tstrand] == "+"
                    qpos = chain[chainidxcur, :qstart] + BP - chain[chainidxcur, :tstart]
                    push!(newpos[idx],
                          (CHR = chain[chainidxcur, :qseq],
                           BP = chain[chainidxcur, :qstrand] == "+" ? qpos : chain[chainidxcur, :qsize] - qpos + 1,
                           strand = chain[chainidxcur, :qstrand],
                           score = chain[chainidxcur, :score]))
                elseif BP < chain[chainidxcur, :tstart]
                    break
                end
            end
            chainidxcur = chainidxcur + 1
        end
    end
    if multiplematches == :error && any(length.(newpos) .> 1)
         error("Error: multiple matches with given positions: $(chrbp[findall(x -> length(x) > 1, newpos), :])")
     elseif multiplematches == :warning && any(length.(newpos) .> 1)
         println("Warning: multiple matches with given positions")
         display(chrbp[findall(x -> length(x) > 1, newpos), :])
     end
     # Return array of lifted over positions (CHR, BP, strand)
     # Same order as arguments
     newpos
end

function findnewcoord(
        CHR::AbstractString, BP::Int64, chain::AbstractDataFrame;
        # Behavior when multiple positions map after liftover
        multiplematches::Symbol = :error # :error, :warning, :silent
    )
    findnewcoord([CHR], Int64[BP], chain; multiplematches = multiplematches)[1]
end

function reversecomplement(seq)
    reversecomplementarr = fill('N', length(seq))
    for i in eachindex(seq)
        if seq[i] == 'A'
            reversecomplementarr[length(seq) - i + 1] = 'T'
        elseif seq[i] == 'T'
            reversecomplementarr[length(seq) - i + 1] = 'A'
        elseif seq[i] == 'C'
            reversecomplementarr[length(seq) - i + 1] = 'G'
        elseif seq[i] == 'G'
            reversecomplementarr[length(seq) - i + 1] = 'C'
        end
    end
    return join(reversecomplementarr, "")
end

function liftoverindel(
        row;
        # pickreference will adjust how reference allele is chosen
        # :first will use the first allele as the reference allele and throw an error 
        # if it does not match
        # :longest will find the longest match and use that as the reference allele
        # :firstlongest will use first allele as the reference allele if it matches, 
        # and if it does not match the reference sequence, use the longest
        pickreference,
        # indelref describes how indel alleles are coded
        # :start assumes that all alleles share the same first nucleotide from the 
        # reference sequence, e.g. A1 = ATCG, A2 = A
        # :startend assumes that all alleles share the same first and last nucleotide 
        # from the reference sequence, e.g. A1 = ATCGA, A2 = AA
        indelref,
        # extendambiguous will extend alleles until they are not substrings of one 
        # another
        # this is what the score plugin does, but if the target and query reference 
        # differ, this will possibly cause the alleles to reflect changes in the reference 
        # as well
        extendambiguous::Bool
    )
    # Make sure that the alleles start and end with the reference sequence
    referenceseq = ""
    alleles = [row.A1, row.A2]
    try
        referenceseq = FASTA.extract(targetfa, row.CHR, pos:(pos + maximum(length.(alleles))))
    catch err
        println(err)
        return nothing
    end
    # Determine which allele to treat as the reference allele
    if pickreference == :first
        if occursin(Regex("^$(alleles[1])"), referenceseq)
            refallele = alleles[1]
        else
            throw("Reference allele does not match reference sequence; pickreference = :first")
        end
    elseif pickreference == :longest
        allelematches =
        [let
             regexmatch = match(Regex("^$(allele)"), referenceseq)
             if !isnothing(regexmatch)
                 regexmatch.match
             else
                 ""
             end
         end
        for allele in alleles]
        refallele = ""
        for allele in allelematches
            if length(allele) > length(refallele)
                refallele = allele
            end
        end
        refallele = string(refallele)
    elseif pickreference == :firstlongest
        if occursin(Regex("^$(alleles[1])"), referenceseq)
            refallele = alleles[1]
        else
        allelematches =
        [let
             regexmatch = match(Regex("^$(allele)"), referenceseq)
             if !isnothing(regexmatch)
                 regexmatch.match
             else
                 ""
             end
         end
        for allele in alleles]
        refallele = ""
        for allele in allelematches
            if length(allele) > length(refallele)
                refallele = allele
            end
        end
        refallele = string(refallele)
        end
    else
        throw("pickreference is invalid, must be one of :first, :longest, :firstlongest")
    end
    reflength = length(refallele)
    if indelref == :start
        # Append an additional nucleotide from the reference sequence
        reflength = reflength + 1
        refend = referenceseq[reflength]
        refallele = refallele * refend
        for i in eachindex(alleles)
            alleles[i] = alleles[i] * refend
        end
    elseif indelref == :startend
        # Check if all the alleles also end with the same nucleotide from the reference sequence
        if !all([allele[end] for allele in alleles] .== refallele[end])
            throw("Alleles' last nucleotide do not match reference sequence, indelref = :startend")
        end
    else
        throw("indelref is invalid, must be one of :start, :startend")
    end
    # Extend sequence if there are any ambiguous alleles
    # If an allele is an exact substring of another, it will make identifying the 
    # reference allele after liftover difficult
    if extendambiguous
        for i in eachindex(alleles)
            for j in (i+1):length(alleles)
                if length(alleles[i]) > length(alleles[j])
                    longerallele = alleles[i]
                    shorterallele = alleles[j]
                elseif length(alleles[i]) < length(alleles[j])
                    longerallele = alleles[j]
                    shorterallele = alleles[i]
                else
                    continue
                end
                while (occursin(Regex("^$(shorterallele)"), longerallele) ||
                    occursin(Regex("$(shorterallele)\$"), longerallele))
                    reflength = reflength + 1
                    refend = ""
                    try
                        refend = FASTA.extract(targetfa, chromosome, (pos+reflength-1):(pos+reflength-1))
                    catch err
                        println(err)
                        return nothing
                    end
                    for k in eachindex(alleles)
                        alleles[k] = alleles[k] * refend
                    end
                    longerallele = longerallele * refend
                    shorterallele = shorterallele * refend
                end
            end
        end
    end
    referenceseq = ""
    try
        referenceseq = FASTA.extract(targetfa, chromosome, pos:(pos + reflength - 1))
    catch err
        println(err)
        return nothing
    end
    # Find the new coordinate
    startcoord = findnewcoord(chromosome, pos, echaindf)
    endcoord = findnewcoord(chromosome, pos + reflength - 1, echaindf)
    if (length(startcoord) == 0 ||
        length(endcoord) == 0 ||
        startcoord[1].chromosome != endcoord[1].chromosome ||
        startcoord[1].strand != endcoord[1].strand)
        return nothing
    end
    startcoord = startcoord[1]
    endcoord = endcoord[1]
    # Flip alleles if needed
    if startcoord.strand == "-"
        qpos = (endcoord.size - endcoord.pos + 1,
                startcoord.size - startcoord.pos + 1)
        for i in eachindex(alleles)
            alleles[i] = reversecomplement.(alleles[i])
        end
    elseif startcoord.strand == "+"
        qpos = (startcoord.pos, endcoord.pos)
    end
    qreferencesequence = ""
    try
        qreferencesequence = FASTA.extract(queryfa, startcoord.chromosome, (qpos[1]):(qpos[2]))
    catch err
        println(err)
        return nothing
    end
    refalleleindex = findall(allele -> allele == qreferencesequence, alleles)
    # Truncate and normalize alleles
    finalend = 0
    while (finalend < minimum(length.(alleles)) &&
           all(qreferencesequence[end - finalend] .==
               [allele[end - finalend] for allele in alleles]))
           finalend = finalend + 1
    end
    qreferencesequence = qreferencesequence[begin:(end - finalend)]
    for i in eachindex(alleles)
        alleles[i] = alleles[i][begin:(end - finalend)]
    end
    finalstart = 1
    while (finalstart < minimum(length.(alleles)) &&
           all(qreferencesequence[finalstart + 1] .==
               [allele[finalstart + 1] for allele in alleles]))
        finalstart = finalstart + 1
    end
    qreferencesequence = qreferencesequence[finalstart:end]
    for i in eachindex(alleles)
        alleles[i] = alleles[i][finalstart:end]
    end
    # Swap alleles to match reference if needed
    return (chromosome = startcoord.chromosome,
            position = qpos[1] + finalstart - 1,
            alleles = alleles,
            refallele = (qreferencesequence, refalleleindex))
end

"""
    liftoversumstats!(gwas::AbstractDataFrame, chain::AbstractDataFrame; kwargs)
    liftoversumstats!(gwas::AbstractVector{<:AbstractDataFrame}, chain::AbstractDataFrame; kwargs)

Perform liftover on a gwas, using an expanded chain file DataFrame 
produced by `GeneticsMakie::readchain`. Variants that are unmapped or 
have multiple matches are dropped. Returns a NamedTuple of DataFrames: 
`unmapped` for unmapped variants (coordinates still on original build), 
`multiple` for variants with multiple matches (coordinates on the 
target build, accompanied by the chain score).
# Arguments
- `multiplematches::Symbol = :error`: Behavior when multiple positions map 
   after liftover. One of :error, :warning, :silent
"""
function liftoversumstats!(
        gwas::AbstractDataFrame,
        chain::AbstractDataFrame;
        # targetfa::FASTX.FASTA.Reader,
        # queryfa::FASTX.FASTA.Reader;
        multiplematches::Symbol = :error,
        # pickreference will adjust how reference allele is chosen
        # :first will use the first allele as the reference allele and throw an error 
        # if it does not match
        # :longest will find the longest match and use that as the reference allele
        # :firstlongest will use first allele as the reference allele if it matches, 
        # and if it does not match the reference sequence, use the longest
        pickreference = :firstlongest,
        # indelref describes how indel alleles are coded
        # :start assumes that all alleles share the same first nucleotide from the 
        # reference sequence, e.g. A1 = ATCG, A2 = A
        # :startend assumes that all alleles share the same first and last nucleotide 
        # from the reference sequence, e.g. A1 = ATCGA, A2 = AA
        indelref = :start,
        # extendambiguous will extend alleles until they are not substrings of one 
        # another
        # this is what the score plugin does, but if the target and query reference 
        # differ, this will possibly cause the alleles to reflect changes in the reference 
        # as well
        extendambiguous::Bool = true
    )
    notlifted = Int64[]
    multiple = Int64[]
    multiplegwas = empty(gwas)
    multiplegwas.score = Int64[]
    newcoords = findnewcoord(gwas[:, :CHR], gwas[:, :BP], chain;
                             multiplematches = multiplematches)
    for i in eachindex(eachrow(gwas))
        if length(newcoords[i]) == 0
            push!(notlifted, i)
            continue
        end
        if length(newcoords[i]) > 1
            additionalmatches = DataFrame(fill(gwas[i, :], length(newcoords[i])))
            for j in eachindex(newcoords[i])
                if newcoords[i][1].strand == "+"
                    additionalmatches[j, :CHR] = newcoords[i][j].CHR
                    additionalmatches[j, :BP] = newcoords[i][j].BP
                elseif newcoords[i][1].strand == "-"
                    if length(gwas[i, :A1]) == 1 || length(gwas[i, :A2]) == 1
                        # SNP
                        gwas[i, :A1] = reversecomplement(gwas[i, :A1])
                        gwas[i, :A2] = reversecomplement(gwas[i, :A2])
                    else
                        # TODO indel
                        push!(notlifted, i)
                        continue
                    end
                end
            end
            additionalmatches[!, :score] = [newcoords[i][j].score
                                            for j in eachindex(newcoords[i])]
            append!(multiplegwas, additionalmatches)
            push!(multiple, i)
            continue
        end
        if newcoords[i][1].strand == "+"
            gwas[i, :CHR] = newcoords[i][1].CHR
            gwas[i, :BP] = newcoords[i][1].BP
        elseif newcoords[i][1].strand == "-"
            if length(gwas[i, :A1]) == 1 || length(gwas[i, :A2]) == 1
                # SNP
                gwas[i, :A1] = reversecomplement(gwas[i, :A1])
                gwas[i, :A2] = reversecomplement(gwas[i, :A2])
            else
                # TODO indel
                push!(notlifted, i)
            end
        end
    end
    unmappedgwas = gwas[notlifted, :]
    delete!(gwas, sort(unique(vcat(notlifted, multiple))))
    (unmapped = unmappedgwas::DataFrame, multiple = multiplegwas::DataFrame)
end

function liftoversumstats!(
        gwas::AbstractVector{<:AbstractDataFrame},
        chain::AbstractDataFrame;
        # targetfa::FASTX.FASTA.Reader,
        # queryfa::FASTX.FASTA.Reader;
        multiplematches::Symbol = :error,
        pickreference = :firstlongest,
        indelref = :start,
        extendambiguous::Bool = true
    )
    arr = liftoversumstats!.(gwas, Ref(chain);
                             multiplematches = multiplematches, pickreference = pickreference, 
                             indelref = indelref, extendambiguous = extendambiguous)
#     arr = liftoversumstats!.(gwas, Ref(chain), Ref(targetfa), Ref(queryfa);
#                              multiplematches = multiplematches, pickreference = pickreference, 
#                              indelref = indelref, extendambiguous = extendambiguous)
    (unmapped = [el.unmapped for el in arr], arr = [el.multiple for el in arr])
end

