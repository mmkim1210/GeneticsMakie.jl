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

function liftoverunmapped(n)
    positions = [Vector{NamedTuple{(:CHR, :BP, :strand, :score), Tuple{String, Int64, String, Int64}}}(undef, 0)
                  for _ in 1:n]
    alleles = [Vector{Vector{NamedTuple{(:alleles, :reference), Tuple{Vector{String}, Int64}}}}(undef, 0)
               for _ in 1:n]
    NamedTuple{(:positions, :alleles), Tuple{Vector{Vector{NamedTuple{(:CHR, :BP, :strand, :score), Tuple{String, Int64, String, Int64}}}}, Vector{Vector{NamedTuple{(:alleles, :reference), Tuple{Vector{String}, Int64}}}}}}((positions = positions, alleles = alleles))
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
    NamedTuple{(:positions, :alleles), Tuple{Vector{Vector{NamedTuple{(:CHR, :BP, :strand, :score), Tuple{String, Int64, String, Int64}}}}, Vector{Vector{NamedTuple{(:alleles, :reference), Tuple{Vector{String}, Int64}}}}}}((positions = snpnewpos, alleles = snpnewalleles))
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
    if whichreference == :first || whichreference == :first_warn || whichreference == :first_silent
        for i in eachindex(referencealleles)
            if !referencealleles[i][1]
                if whichreference == :first
                    throw("$(CHR[i]):$(BP[i]) reference allele does not match reference sequence")
                elseif whichreference == :first_warn
                    println("$(CHR[i]):$(BP[i]) reference allele does not match reference sequence")
                    referencelengths[i] = 0
                    continue
                elseif whichreference == :first_silent
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
    # Breaks if the reference sequence length is longer than 1000, which is a bit of 
    # an arbitrary limit, but chances are that something has gone wrong
    if extendambiguous
        for i in eachindex(indelalleles)
            if (referencelengths[i] == 0 || referencelengths[i] > 1000 ||
                length(unique(indelalleles[i])) != length(indelalleles[i]))
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
                        if referencelengths[i] > 1000
                            break
                        end
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
                  if newreferencealleles[i] == 0
                      nothing
                  else
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
                      if startidx == 0
                          newreferencealleles[i] = 0
                          nothing
                      else
                          for j in eachindex(liftedoveralleles[i])
                              liftedoveralleles[i][j] = liftedoveralleles[i][j][startidx:end]
                          end
                          qreferencesequence = qreferencesequence[startidx:end]
                          (CHR = bounds[1].CHR, BP = bounds[1].BP + startidx - 1, strand = bounds[1].strand, 
                           score = bounds[1].score)
                      end
                  end
              end
              for (i, bounds) in enumerate(possiblebounds)]
             badidxs = findall(idx -> idx == 0, newreferencealleles)
             deleteat!(indelnewpos, badidxs)
             deleteat!(liftedoveralleles, badidxs)
             deleteat!(newreferencealleles, badidxs)
             (positions = indelnewpos,
              alleles = [(alleles = liftedoveralleles[i],
                          reference = newreferencealleles[i])
                         for i in eachindex(liftedoveralleles)])
         end
     end
     for i in eachindex(CHR)] |>
    liftedindels -> map(fieldname -> getfield.(liftedindels, fieldname), [:positions, :alleles])
    NamedTuple{(:positions, :alleles), Tuple{Vector{Vector{NamedTuple{(:CHR, :BP, :strand, :score), Tuple{String, Int64, String, Int64}}}}, Vector{Vector{NamedTuple{(:alleles, :reference), Tuple{Vector{String}, Int64}}}}}}((positions = liftedpositions, alleles = liftedalleles))
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
    badidxs = findall([all(.!occursin.(r"^[ATCG]*$", posalleles)) for posalleles in alleles])
    goodidxs = setdiff(eachindex(alleles), badidxs)
    snpidxs = findall([all(length.(alleles[idx]) .== 1) for idx in goodidxs])
    indelidxs = setdiff(goodidxs, snpidxs)
    # Output unmapped positions with malformed alleles
    badnewpos, badnewalleles =
    liftoverunmapped(length(badidxs))
    # Liftover SNPs
    snpnewpos, snpnewalleles =
    liftoversnps(CHR[snpidxs], BP[snpidxs], alleles[snpidxs], targetfa, queryfa, chain;
                 multiplematches = multiplematches)
    # Liftover indels
    indelnewpos, indelnewalleles =
    liftoverindels(CHR[indelidxs], BP[indelidxs], alleles[indelidxs], targetfa, queryfa, 
                 chain;
                 multiplematches = multiplematches, whichreference = whichreference, 
                 extendambiguous = extendambiguous)

    positions = Vector{Vector{NamedTuple{(:CHR, :BP, :strand, :score), Tuple{String, Int64, String, Int64}}}}(undef, length(CHR))
    alleles = Vector{Vector{NamedTuple{(:alleles, :reference), Tuple{Vector{String}, Int64}}}}(undef, length(CHR))
    for (i, idx) in enumerate(badidxs)
        positions[idx] = badnewpos[i]
        alleles[idx] = badnewalleles[i]
    end
    for (i, idx) in enumerate(snpidxs)
        positions[idx] = snpnewpos[i]
        alleles[idx] = snpnewalleles[i]
    end
    for (i, idx) in enumerate(indelidxs)
        positions[idx] = indelnewpos[i]
        alleles[idx] = indelnewalleles[i]
    end
    NamedTuple{(:positions, :alleles), Tuple{Vector{Vector{NamedTuple{(:CHR, :BP, :strand, :score), Tuple{String, Int64, String, Int64}}}}, Vector{Vector{NamedTuple{(:alleles, :reference), Tuple{Vector{String}, Int64}}}}}}((positions = positions, alleles = alleles))
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
   :first_silent will use the first allele as the reference allele and silently mark 
   the position unmapped if it does not match
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
    @assert whichreference ∈ [:first, :first_warn, :first_silent, :longest]
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

