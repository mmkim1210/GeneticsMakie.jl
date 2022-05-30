"""
    parsegtf!(gencode::DataFrame)

Parse `gencode` by extracting `gene_id`, `gene_name`, `gene_type`, `transcript_id`,
`transcript_support_level` information from the `info` column.
"""
function parsegtf!(gencode::DataFrame)
    if !("info" in names(gencode))
        error("Cannot find :info column in the annotation.")
    end
    for col in [:gene_id, :gene_name, :gene_type]
        gencode[:, col] =
            [getindex(m.captures, 1) for m in match.(Regex("$(col) \"(.*?)\";"), gencode.info)]
    end
    for col in [:transcript_id, :transcript_support_level]
        storage = Vector(undef, size(gencode, 1))
        for (n, m) in enumerate(match.(Regex("$(col) \"(.*?)\";"), gencode.info))
            if isnothing(m) && col == :transcript_id
                storage[n] = gencode.gene_id[n]
            elseif isnothing(m) && col == :transcript_support_level
                storage[n] = missing
            else
                storage[n] = getindex(m.captures, 1)
            end
        end
        gencode[:, col] = storage
    end
    gencode.gene_id = [getindex(i, 1) for i in split.(gencode.gene_id, ".")]
    gencode.transcript_id = [getindex(i, 1) for i in split.(gencode.transcript_id, ".")]
    gencode.seqnames = replace.(gencode.seqnames, "chr" => "")
    return
end

"""
    findgene(gene::AbstractString, gencode::DataFrame)

Find chromosome, gene start, and gene stop sites for the `gene` of interest.
"""
function findgene(gene::AbstractString, gencode::DataFrame)
    if startswith(gene, "ENSG")
        ind = findfirst(isequal(gene), gencode.gene_id)
    else
        ind = findfirst(isequal(gene), gencode.gene_name)
    end
    if isnothing(ind)
        error("Cannot find $(gene) in the annotation.")
    else
        return gencode.seqnames[ind], gencode.start[ind], gencode[ind, :end]
    end
end