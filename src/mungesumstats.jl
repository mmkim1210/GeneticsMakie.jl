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
                    # Writing code to handle "-" strand could be possible but probably 
                    # wouldn't make sense to do; a proper chain file should have all 
                    # the blocks for the build we're lifting from on the + strand
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

# Wrapper for FASTA.extract that will catch errors and print warnings + other convenience 
# functions
function FASTAextract_warn(fasta, CHR, BPrange)
    chromosomenames = keys(fasta.index.names)
    try
        if CHR ∈ chromosomenames
            FASTA.extract(fasta, CHR, BPrange)
        elseif replace(CHR, r"^chr" => "") ∈ chromosomenames
            # Handle differences in chromosome naming
            FASTA.extract(fasta, replace(CHR, r"^chr" => ""), BPrange)
        elseif "chr" * CHR ∈ chromosomenames
            FASTA.extract(fasta, "chr" * CHR, BPrange)
        else
            nothing
        end
    catch err
        println(err)
        nothing
    end
end

function snpreference(
        CHR::Vector{<:AbstractString}, BP::Vector{Int64}, alleles::Vector{<:Vector{<:AbstractString}},
        queryfa)
    @assert length(CHR) == length(BP) == length(alleles)
    [let
         refseq = FASTAextract_warn(queryfa, CHR[i], BP[i]:BP[i])
         findfirst(allele -> allele == refseq, alleles[i]) |>
         refidx -> isnothing(refidx) ? 0 : refidx
     end
     for i in eachindex(CHR)]
end

function snpreference(
        CHR::AbstractString, BP::Int64, alleles::Vector{<:AbstractString},
        queryfa)
    snpreference([CHR], [BP], [alleles], queryfa)[1]
end

function liftoversnps(
        CHR::Vector{<:AbstractString}, BP::Vector{Int64}, alleles::Vector{<:Vector{<:AbstractString}},
        targetfa, queryfa, chain::AbstractDataFrame;
        # Behavior when multiple positions map after liftover
        # :error, :warning, :silent
        multiplematches::Symbol = :error,
    )
    snpnewpos = findnewcoord(CHR, BP, chain;
                             multiplematches = multiplematches)
    snpnewalleles=
    [[let
          if snpnewpos[i][j].strand == "+"
              snpnewalleles = deepcopy(alleles[i])
          elseif snpnewpos[i][j].strand == "-"
              snpnewalleles = reversecomplement.(alleles[i])
          end
          snpnewreference =
          snpreference(snpnewpos[i][j].CHR, snpnewpos[i][j].BP, snpnewalleles, 
                       queryfa)
          (alleles = snpnewalleles, reference = snpnewreference)
      end
      for j in eachindex(snpnewpos[i])]
     for i in eachindex(snpnewpos)]
    (positions = snpnewpos, alleles = snpnewalleles)
end

function liftoverindels(
        CHR::Vector{<:AbstractString}, BP::Vector{Int64}, alleles::Vector{<:Vector{<:AbstractString}},
        targetfa, queryfa, chain::AbstractDataFrame;
        # Behavior when multiple positions map after liftover
        # :error, :warning, :silent
        multiplematches::Symbol = :error,
        whichreference::Symbol = :first,
        indelreference::Symbol = :start,
        extendambiguous::Bool = true
    )
    # Liftover indels
    # Copy the alleles so that they are editable
    indelalleles = deepcopy(alleles)
    # Find reference allele
    referenceseqs =
    [FASTAextract_warn(targetfa,
                       CHR[i],
                       BP[i]:(BP[i] + maximum(length.(alleles[i]))))
     for i in eachindex(CHR)]
    referencealleles =
    [isnothing(referenceseqs[i]) ? false : startswith.(referenceseqs[i], indelalleles[i])
     for i in eachindex(CHR)]
    referencelengths = fill(0, length(referencealleles))
    if whichreference == :first || whichreference == :first_warn
        for i in eachindex(referencealleles)
            if !referencealleles[i][1]
                if whichreference == :first
                    throw("$(CHR[i]):$(BP[i]) reference allele does not match reference sequence")
                elseif whichreference == :first_warn
                    println("$(CHR[i]):$(BP[i]) reference allele does not match reference sequence")
                    referencelengths[i] = 0
                    continue
                end
            end
            referencealleles[i][2:end] .= false
            referencelengths[i] = length(indelalleles[i][1])
        end
    elseif whichreference == :longest
        for i in eachindex(referencealleles)
            longestallele = 0
            longestallelelength = 0
            for j in eachindex(referencealleles[i])
                if referencealleles[i][j]
                    curlength = length(indelalleles[i][j])
                    if curlength > longestallelelength
                        longestallele = j
                        longestalellelength = curlength
                    end
                end
            end
            referencealleles[i] .= false
            longestallele != 0 ? referencealleles[i][longestallele] = true : nothing
            referencelengths[i] = longestallelelength
        end
    end
    # Prepare and check the indelalleles array
    if indelreference == :start
        # Append an additional nucleotide from the reference sequence
        for i in eachindex(referencelengths)
            if referencelengths[i] != 0
                referencelengths[i] = referencelengths[i] + 1
            end
        end
        for i in eachindex(indelalleles)
            if referencelengths[i] == 0
                continue
            end
            refend = referenceseqs[i][referencelengths[i]]
            for j in eachindex(indelalleles[i])
                indelalleles[i][j] = indelalleles[i][j] * refend
            end
        end
    elseif indelref == :startend
        # Check if all the alleles also end with the same nucleotide from the reference sequence
        for i in eachindex(indelalleles)
            if referencelengths[i] == 0
                continue
            end
            if (!all([allele[end] for allele in indelalleles[i]] .== referenceseqs[i][referencelengths[i]]) ||
                !all([allele[begin] for allele in indelalleles[i]] .== referenceseqs[i][begin]))
                throw("Alleles' surrounding nucleotides do not match reference sequence")
            end
        end
    end
    # Extend sequence if there are any ambiguous alleles
    # If an allele is an exact substring of another, it will make identifying the 
    # reference allele after liftover difficult
    if extendambiguous
        for i in eachindex(indelalleles)
            if referencelengths[i] == 0
                continue
            end
            for j in eachindex(indelalleles[i])
                for k in (j+1):length(indelalleles[i])
                    if length(indelalleles[i][j]) > length(indelalleles[i][k])
                        longerallele = indelalleles[i][j]
                        shorterallele = indelalleles[i][k]
                    elseif length(indelalleles[i][j]) < length(indelalleles[i][k])
                        shorterallele = indelalleles[i][j]
                        longerallele = indelalleles[i][k]
                    else
                        continue
                    end
                    while (occursin(Regex("^$(shorterallele)"), longerallele) ||
                           occursin(Regex("$(shorterallele)\$"), longerallele))
                        referencelengths[i] = referencelengths[i] + 1
                        refend = ""
                        refend =
                        FASTAextract_warn(targetfa,
                                          CHR[i],
                                          (BP[i]+referencelengths[i]-1):(BP[i]+referencelengths[i]-1))
                        for l in eachindex(indelalleles[i])
                            indelalleles[i][l] = indelalleles[i][l] * refend
                        end
                        longerallele = indelalleles[i][j]
                        shorterallele = indelalleles[i][k]
                    end
                end
            end
        end
    end
    referenceseqs =
    [referencelengths[i] != 0 ? FASTAextract_warn(targetfa,
                                                  CHR[i],
                                                  BP[i]:(BP[i] + referencelengths[i] - 1)) : nothing
     for i in eachindex(CHR)]
    # Find the new coordinate
    startcoords = findnewcoord(CHR, BP, chain;
                               multiplematches = multiplematches)
    endcoords = findnewcoord(CHR, BP .+ referencelengths .- 1, chain;
                             multiplematches = multiplematches)
    liftedpositions, liftedalleles =
    [let
         if referencelengths[i] == 0
             (positions = Vector{NamedTuple{(:CHR, :BP, :strand, :score), Tuple{String, Int64, String, Int64}}}(undef, 0),
              alleles = Vector{Vector{NamedTuple{(:alleles, :reference), Tuple{Vector{String}, Int64}}}}(undef, 0))
         else
             possiblebounds =
             Iterators.product(startcoords[i], endcoords[i]) |>
             iter -> Iterators.flatten((iter,))
             liftedoverbounds =
             Iterators.filter(startend -> (startend[1].CHR == startend[2].CHR &&
                                           abs(startend[1].BP - startend[2].BP) + 1 == referencelengths[i] &&
                                           startend[1].strand == startend[2].strand),
                              possiblebounds)
             liftedoveralleles =
             [startend[1].strand == "+" ? deepcopy(indelalleles[i]) : reversecomplement.(indelalleles[i])
              for startend in possiblebounds]
             newreferencealleles = fill(0, length(liftedoveralleles))
             indelnewpos =
             [let
                  # Left aligning alleles
                  qreferencesequence = ""
                  qreferencesequence = FASTAextract_warn(queryfa,
                                                         bounds[1].CHR,
                                                         (bounds[1].BP):(bounds[2].BP))
                  newreferencealleles[i] =
                  findfirst(seq -> seq == qreferencesequence, liftedoveralleles[i]) |>
                  refidx -> isnothing(refidx) ? 0 : refidx
                  endidx = minimum(length.(vcat(liftedoveralleles[i], qreferencesequence))) - 1
                  negidx = 0
                  while ((negidx < endidx) &&
                         (all(qreferencesequence[end - negidx] .==
                              [allele[end - negidx] for allele in liftedoveralleles[i]])))
                      negidx = negidx + 1
                  end
                  for j in eachindex(liftedoveralleles[i])
                      liftedoveralleles[i][j] = liftedoveralleles[i][j][begin:(end-negidx)]
                  end
                  qreferencesequence = qreferencesequence[begin:(end-negidx)]
                  startidx = 0
                  endidx = minimum(length.(vcat(liftedoveralleles[i], qreferencesequence)))
                  while ((startidx < endidx) &&
                         all(qreferencesequence[begin + startidx] .==
                             [allele[begin + startidx] for allele in liftedoveralleles[i]]))
                      startidx = startidx + 1
                  end
                  for j in eachindex(liftedoveralleles[i])
                      liftedoveralleles[i][j] = liftedoveralleles[i][j][startidx:end]
                  end
                  qreferencesequence = qreferencesequence[startidx:end]
                  (CHR = bounds[1].CHR, BP = bounds[1].BP + startidx - 1, strand = bounds[1].strand, 
                   score = bounds[1].score)
              end
              for (i, bounds) in enumerate(possiblebounds)]
             (positions = indelnewpos,
              alleles = [(alleles = liftedoveralleles[i],
                          reference = newreferencealleles[i])
                         for i in eachindex(liftedoveralleles)])
         end
     end
     for i in eachindex(CHR)] |>
    liftedindels -> map(fieldname -> getfield.(liftedindels, fieldname), [:positions, :alleles])
    (positions = liftedpositions, alleles = liftedalleles)
end

function liftoveralleles(
        CHR::Vector{<:AbstractString}, BP::Vector{Int64}, alleles::Vector{<:Vector{<:AbstractString}},
        targetfa, queryfa, chain::AbstractDataFrame;
        # Behavior when multiple positions map after liftover
        # :error, :warning, :silent
        multiplematches::Symbol = :error,
        whichreference::Symbol = :first,
        indelreference::Symbol = :start,
        extendambiguous::Bool = true
    )
    # Split into SNPs and indels
    snpidxs = findall([all(length.(posalleles) .== 1) for posalleles in alleles])
    indelidxs = setdiff(eachindex(alleles), snpidxs)

    # Liftover SNPs
    snpnewpos, snpnewalleles =
    liftoversnps(CHR[snpidxs], BP[snpidxs], alleles[snpidxs], targetfa, queryfa, chain,
                 multiplematches = multiplematches)

    # Liftover indels
    indelnewpos, indelnewalleles =
    liftoverindels(CHR[indelidxs], BP[indelidxs], alleles[indelidxs], targetfa, queryfa, 
                 chain;
                 multiplematches = multiplematches, whichreference = whichreference, 
                 extendambiguous = extendambiguous)

    positions = Vector{Vector{NamedTuple{(:CHR, :BP, :strand, :score), Tuple{String, Int64, String, Int64}}}}(undef, length(CHR))
    alleles = Vector{Vector{NamedTuple{(:alleles, :reference), Tuple{Vector{String}, Int64}}}}(undef, length(CHR))
    for (i, idx) in enumerate(snpidxs)
        positions[idx] = snpnewpos[i]
        alleles[idx] = snpnewalleles[i]
    end
    for (i, idx) in enumerate(indelidxs)
        positions[idx] = indelnewpos[i]
        alleles[idx] = indelnewalleles[i]
    end
    (positions = positions, alleles = alleles)
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
- `whichreference::Symbol = :first`: For indels, describes how the reference 
   allele is chosen. This is important when the alelles have to be extended using 
   the reference genome, which is relevant when the lifted over segment is on the 
   opposite strand.  
   :first will use the first allele as the reference allele and throw an error if 
   it does not match  
   :first_warn will use the first allele as the reference allele and warn the user 
   if it does not match; the position will be unmapped
   :longest will find the allele with the longest match and use that as the reference 
   allele  
- `indelreference::Symbol = :start`: For indels, describes how indel alleles are coded.
   :start assumes that all alleles share the same first nucleotide from the 
   reference sequence, e.g. A1 = ATCG, A2 = A  
   :startend assumes that all alleles share the same first and last nucleotide 
   from the reference sequence, e.g. A1 = ATCGA, A2 = AA
- `extendambiguous::Bool = true`: Extend alleles until they are not substrings of 
   one another.
- `referenceorder::Bool = true`: Reorder alleles so that referenece allele is A1.
"""
function liftoversumstats!(
        gwas::AbstractDataFrame,
        targetfa::FASTA.Reader,
        queryfa::FASTA.Reader,
        chain::AbstractDataFrame;
        multiplematches::Symbol = :error,
        whichreference::Symbol = :first,
        indelreference = :start,
        extendambiguous::Bool = true,
        referenceorder::Bool = true
    )
    @assert multiplematches ∈ [:error, :warning, :silent]
    @assert whichreference ∈ [:first, :first_warn, :longest]
    @assert indelreference ∈ [:start, :startend]
    # Indices for rows that did not properly liftover
    notlifted = Int64[]
    multiple = Int64[]
    # DataFrame for rows that mapped to multiple positions
    multiplegwas = empty(gwas)
    multiplegwas.score = Int64[]

    newpositions, newalleles =
    liftoveralleles(gwas[:, :CHR], gwas[:, :BP], [[row.A1, row.A2] for row in eachrow(gwas)],
                    targetfa, queryfa, chain;
                    multiplematches = multiplematches, whichreference = whichreference, 
                    indelreference = indelreference, extendambiguous = extendambiguous)
    for i in eachindex(eachrow(gwas))
        if length(newpositions[i]) == 0
            push!(notlifted, i)
            continue
        elseif length(newpositions[i]) == 1
            gwas[i, :CHR] = newpositions[i][1].CHR
            gwas[i, :BP] = newpositions[i][1].BP
            if referenceorder
                if newalleles[i][1].reference == 1
                    gwas[i, :A1] = newalleles[i][1].alleles[1]
                    gwas[i, :A2] = newalleles[i][1].alleles[2]
                elseif newalleles[i][1].reference == 2
                    gwas[i, :A1] = newalleles[i][1].alleles[2]
                    gwas[i, :A2] = newalleles[i][1].alleles[1]
                    if :Z ∈ propertynames(gwas)
                        gwas[i, :Z] = -gwas[i, :Z]
                    end
                    if :BETA ∈ propertynames(gwas)
                        gwas[i, :BETA] = -gwas[i, :BETA]
                    end
                    if :OR ∈ propertynames(gwas)
                        gwas[i, :OR] = 1 / gwas[i, :OR]
                    end
                else
                    push!(notlifted, i)
                    continue
                end
            else
                gwas[i, :A1] = newalleles[i][1].alleles[1]
                gwas[i, :A2] = newalleles[i][1].alleles[2]
            end
        elseif length(newcoords[i]) > 1
            additionalmatches = DataFrame(fill(gwas[i, :], length(newpositions[i])))
            referencemismatch = []
            for j in eachindex(newpositions[i])
                additionalmatches[j, :CHR] = newpositions[i][j].CHR
                additionalmatches[j, :BP] = newpositions[i][j].BP
                if referenceorder
                    if newalleles[i][j].reference == 1
                        additionalmatches[j, :A1] = newalleles[i][j].alleles[1]
                        additionalmatches[j, :A2] = newalleles[i][j].alleles[2]
                    elseif newalleles[i][j].reference == 2
                        additionalmatches[j, :A1] = newalleles[i][j].alleles[2]
                        additionalmatches[j, :A2] = newalleles[i][j].alleles[1]
                        if :Z ∈ propertynames(additionalmatches)
                            additionalmatches[j, :Z] = -additionalmatches[j, :Z]
                        end
                        if :BETA ∈ propertynames(additionalmatches)
                            additionalmatches[j, :BETA] = -additionalmatches[j, :BETA]
                        end
                        if :OR ∈ propertynames(additionalmatches)
                            additionalmatches[j, :OR] = 1 / additionalmatches[j, :OR]
                        end
                    else
                        push!(referencemismatch, j)
                        continue
                    end
                else
                    additionalmatches[j, :A1] = newalleles[i][j].alleles[1]
                    additionalmatches[j, :A2] = newalleles[i][j].alleles[2]
                end
            end
            additionalmatches[!, :score] = [newpositions[i][j].score
                                            for j in eachindex(newcoords[i])]
            if referenceorder
                deleteat!(additionalmatches, referencemismatch)
            end
            append!(multiplegwas, additionalmatches)
            push!(multiple, i)
            continue
        end
    end
    unmappedgwas = gwas[notlifted, :]
    delete!(gwas, sort(unique(vcat(notlifted, multiple))))
    (unmapped = unmappedgwas::DataFrame, multiple = multiplegwas::DataFrame)
end

function liftoversumstats!(
        gwas::AbstractVector{<:AbstractDataFrame},
        targetfa::FASTA.Reader,
        queryfa::FASTA.Reader,
        chain::AbstractDataFrame;
        multiplematches::Symbol = :error,
        whichreference = :first,
        indelreference = :start,
        extendambiguous::Bool = true,
        referenceorder::Bool = true
    )
    arr = liftoversumstats!.(gwas, Ref(targetfa), Ref(queryfa), Ref(chain);
                             multiplematches = multiplematches, whichreference = whichreference, 
                             indelreference = indelreference, extendambiguous = extendambiguous, 
                             referenceorder = referenceorder)
    (unmapped = [el.unmapped for el in arr], arr = [el.multiple for el in arr])
end

