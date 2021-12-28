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

function findgwasloci(gwas::DataFrame; p = 5e-8)
    loci = DataFrame(CHR = String[], BP = Int[], P = Float64[])
    df = gwas[findall(x -> x < p, gwas.P), ["CHR", "BP", "P"]]
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

function findsnps(CHR₁::AbstractVector, BP₁::AbstractVector, A1₁::AbstractVector, A2₁::AbstractVector, 
    CHR₂::AbstractVector, BP₂::AbstractVector, A1₂::AbstractVector, A2₂::AbstractVector)
    
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
end

findsnps(gwas::DataFrame, ref::SnpData) = findsnps(gwas.CHR, gwas.BP, gwas.A1, gwas.A2, 
    ref.snp_info.chromosome, ref.snp_info.position, ref.snp_info.allele1, ref.snp_info.allele2)

findsnps(ref::SnpData, gwas::DataFrame) = 
    findsnps(ref.snp_info.chromosome, ref.snp_info.position, ref.snp_info.allele1, ref.snp_info.allele2,
        gwas.CHR, gwas.BP, gwas.A1, gwas.A2)

findsnps(gwas₁::DataFrame, gwas₂::DataFrame) = findsnps(gwas₁.CHR, gwas₁.BP, gwas₁.A1, gwas₁.A2, gwas₂.CHR, gwas₂.BP, gwas₂.A1, gwas₂.A2)

function findmissing(ind::Matrix{Union{Missing, Int64}})
    storage = Vector{Union{Missing, Int}}(undef, size(ind, 1))
    for i in eachindex(storage)
        if sum(ismissing.(ind[i, :])) == 2
            storage[i] = missing
        elseif ismissing(ind[i, 1])
            storage[i] = ind[i, 2]
        else
            storage[i] = ind[i, 1]
        end
    end
    return storage
end

gwas = Dict(
    "scz" => (url = "https://figshare.com/ndownloader/files/28169757",
        PMID = "", title = "Schizophrenia (PGC3)", file = "PGC3_SCZ_wave3_public.v2.tsv.gz"),
    "bd" => (url = "https://figshare.com/ndownloader/files/26603681",
        PMID = "34002096", title = "Bipolar (Mullins et al. 2021)", file = "pgc-bip2021-all.vcf.tsv.gz"),
    "asd" => (url = "https://figshare.com/ndownloader/files/28169292",
        PMID = "30804558", title = "Autism (Grove et al. 2019)", file = "iPSYCH-PGC_ASD_Nov2017.gz"),
    "scz vs bd" => (url = "https://figshare.com/ndownloader/files/28169358",
        PMID = "29906448", title = "SCZ vs BD (Ruderfer et al. 2018)", file = "SCZvsBD.sumstats.gz"),
    "adhd" => (url = "https://figshare.com/ndownloader/files/28169253",
        PMID = "30478444", title = "ADHD (Demontis et al. 2019)", file = "daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz"),
    "anorexia" => (url = "https://figshare.com/ndownloader/files/28169271",
        PMID = "31308545", title = "Anorexia (Watson et al. 2019)", file = "pgcAN2.2019-07.vcf.tsv.gz"),
    "amd" => (url = "http://csg.sph.umich.edu/abecasis/public/amd2015/Fritsche_2015_AdvancedAMD.txt.gz",
        PMID = "26691988", title = "Age-related macular degeneration (Fritsche et al. 2016)", file = "Fritsche_2015_AdvancedAMD.txt.gz"),
    "menarche" => (url = "https://www.reprogen.org/Menarche_1KG_NatGen2017_WebsiteUpload.zip",
        PMID = "28436984", title = "Age at menarche (Day et al. 2017)", file = "Menarche_1KG_NatGen2017_WebsiteUpload.txt"),
    "menopause" => (url = "https://www.reprogen.org/reprogen_ANM_201K_170621.txt.gz",
        PMID = "34349265", title = "Age at menoapuse (Ruth et al. 2021)", file = "reprogen_ANM_201K_170621.txt.gz"),
    "alz" => (url = "https://ctg.cncr.nl/documents/p1651/AD_sumstats_Jansenetal_2019sept.txt.gz",
        PMID = "30617256", title = "Alzheimer disease (Jansen et al. 2019)", file = "AD_sumstats_Jansenetal_2019sept.txt.gz"),
    "cross" => (url = "https://figshare.com/ndownloader/files/28169382",
        PMID = "31835028", title = "Cross disorder (PGC 2019)", file = "pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt.gz"),
    "mdd" => (url = "https://datashare.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt?sequence=3&isAllowed=y",
        PMID = "30718901", title = "Major depression (Howard et al. 2019)", file = "PGC_UKB_depression_genome-wide.txt"),
    "lupus" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003156/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz",
        PMID = "26502338", title = "Lupus (Bentham et al. 2015)", file = "bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz"),
    "height" => (url = "https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz",
        PMID = "30124842", title = "Height (Yengo et al. 2018)", file = "Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz"),
    "weight" => (url = "https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz",
        PMID = "30124842", title = "Weight (Yengo et al. 2018)", file = "Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz"),
    "left-hand" => (url = "https://evansgroup.di.uq.edu.au/GWAS_RESULTS/HANDEDNESS/LeftHandedness_MetaAnalysis_UKBB_IHC.txt.gz",
        PMID = "32989287", title = "Left-handedness (Cuellar-Partida et al. 2020)", file = "LeftHandedness_MetaAnalysis_UKBB_IHC.txt.gz"),
    "tsh" => (url = "http://csg.sph.umich.edu/willer/public/TSH2020/2020_Zhou_et_al_TSH_meta_HUNT-MGI-ThyroidOmics.txt",
        PMID = "32769997", title = "TSH level (Zhou et al. 2020)", file = "2020_Zhou_et_al_TSH_meta_HUNT-MGI-ThyroidOmics.txt"),
    "ra" => (url = "http://plaza.umin.ac.jp/~yokada/datasource/files/GWASMetaResults/RA_GWASmeta_European_v2.txt.gz",
        PMID = "24390342", title = "Rheumatoid arthritis (Okada et al. 2014)", file = "RA_GWASmeta_European_v2.txt.gz"),
    "als" => (url = "https://surfdrive.surf.nl/files/index.php/s/E5RetKw10hC3jXy/download?path=%2F&files=ALS_GWAS_RVB_SMR_sumstats.tar.gz",
        PMID = "34873335", title = "Amyotrophic lateral sclerosis (van Rheenen et al. 2021)", file = "ALS_sumstats_EUR_only.txt.gz"),
    "parkinson" => (url = "https://drive.google.com/file/d/1FZ9UL99LAqyWnyNBxxlx6qOUlfAnublN/view?usp=sharing",
        PMID = "31701892", title = "Parkinson's disease (Nalls et al. 2019)", file = "nallsEtAl2019_excluding23andMe_allVariants.tab"),
    "cp" => (url = "https://thessgac.com/papers/", PMID = "30038396", title = "Cognitive performance (Lee et al. 2018)", file = "GWAS_CP_all.txt"),
    "ea" => (url = "https://thessgac.com/papers/", PMID = "30038396", title = "Educational attainment (Lee et al. 2018)", file = "GWAS_EA_excl23andMe.txt"),
    "birth" => (url = "https://thessgac.com/papers/", PMID = "27798627", title = "Age at first birth (Barban et al. 2016)", file = "AgeFirstBirth_Pooled.txt"),
    "children" => (url = "https://thessgac.com/papers/", PMID = "27798627", title = "Number of children (Barban et al. 2016)", file = "NumberChildrenEverBorn_Pooled.txt"),
    "risk" => (url = "https://thessgac.com/papers/", PMID = "30643258", title = "Risk tolerance (Linner et al. 2019)", file = "RISK_GWAS_MA_UKB+replication.txt"),
    "neuroticism" => (url = "https://ctg.cncr.nl/documents/p1651/sumstats_neuroticism_ctg_format.txt.gz",
        PMID = "29942085", title = "Neuroticism (Nagel et al. 2018)", file = "sumstats_neuroticism_ctg_format.txt.gz"),
    "asthma" => (url = "https://genepi.qimr.edu.au/staff/manuelf/gwas_results/SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz",
        PMID = "29083406", title = "Asthma (Ferreira et al. 2017)", file = "SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz"),
    "cad" => (url = "http://www.cardiogramplusc4d.org/data-downloads/",
        PMID = "26343387", title = "CAD (Nikpay et al. 2015)", file = "cad.add.160614.website.txt"),
    "afib" => (url = "https://personal.broadinstitute.org/ryank/AF_HRC_GWAS_ALLv11.zip",
        PMID = "29892015", title = "Atrial fibrillation (Roselli et al. 2018)", file = "AF_HRC_GWAS_ALLv11.txt"),
    "stroke" => (url = "https://www.kp4cd.org/datasets/stroke",
        PMID = "29531354", title = "Stroke (Malik et al. 2018)", file = "metastroke.all.chr.bp"),
    "ibd" => (url = "https://www.ibdgenetics.org/downloads.html",
        PMID = "26192919", title = "Inflammatory bowel disease (Liu et al. 2015)", file = "EUR.IBD.gwas_info03_filtered.assoc"),
    "ckd" => (url = "http://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
        PMID = "31152163", title = "Chronic kidney disease (Wuttke et al. 2019)", file = "CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz"),
    "t2d" => (url = "https://diagram-consortium.org/downloads.html",
        PMID = "30297969", title = "Type II diabetes (Mahajan et al. 2018)", file = "Mahajan.NatGenet2018b.T2D.European.txt"),
    "smoking" => (url = "https://conservancy.umn.edu/bitstream/handle/11299/201564/CigarettesPerDay.txt.gz",
        PMID = "30643251", title = "Cigarettes per day (Liu et al. 2019)", file = "CigarettesPerDay.txt.gz"),
    "alcohol" => (url = "https://conservancy.umn.edu/bitstream/handle/11299/201564/DrinksPerWeek.txt.gz",
        PMID = "30643251", title = "Drinks per week (Liu et al. 2019)", file = "DrinksPerWeek.txt.gz"),
    "intelligence" => (url = "https://ctg.cncr.nl/documents/p1651/SavageJansen_IntMeta_sumstats.zip",
        PMID = "29942086", title = "Intelligence (Savage et al. 2018)", file = "SavageJansen_2018_intelligence_metaanalysis.txt"),
    "insomnia" => (url = "https://ctg.cncr.nl/documents/p1651/Insomnia_sumstats_Jansenetal.txt.gz",
        PMID = "30804565", title = "Insomnia (Jansen et al. 2019)", file = "Insomnia_sumstats_Jansenetal.txt.gz"),
    "icv" => (url = "https://ctg.cncr.nl/documents/p1651/meta_analysis_BV_Jansenetal_2020.sumstats.txt.gz",
        PMID = "33154357", title = "Intracranial volume (Jansen et al. 2020)", file = "meta_analysis_BV_Jansenetal_2020.sumstats.txt.gz"),
    "covid" => (url = "https://storage.googleapis.com/covid19-hg-public/20210415/results/20210607/COVID19_HGI_A2_ALL_leave_23andme_20210607.b37.txt.gz",
        PMID = "34237774", title = "COVID-19 (2021)", file = "COVID19_HGI_A2_ALL_leave_23andme_20210607.b37.txt.gz")
)
# add glaucoma, copd, prostate/breast/colorectal/lung cancers, celiac disease, ptsd, male-pattern baldness,
# multiple sclerosis, hypertension, neuroimaging measures, beat synchronization

function downloadgwas(path::AbstractString; pheno::Union{AbstractVector, AbstractString, Nothing} = nothing)
    if isnothing(pheno)
        for key in keys(gwas)
            key in ["als", "menarche", "parkinson", "intelligence", "t2d", "ibd", "stroke", "afib", "cad", "asthma",
                "cp", "ea", "birth", "children", "risk"] ? continue : nothing
            @info("Downloading summary statistics for $(key)")
            url = gwas[key].url
            file = gwas[key].file
            isdir(path) || mkdir(path)
            isfile("$(path)/$(file)") || download(url, "$(path)/$(file)")
        end
    elseif isa(pheno, AbstractVector)
        for key in pheno
            @info("Downloading summary statistics for $(key)")
            url = gwas[key].url
            file = gwas[key].file
            isdir(path) || mkdir(path)
            isfile("$(path)/$(file)") || download(url, "$(path)/$(file)")
        end
    else
        @info("Downloading summary statistics for $(pheno)")
        url = gwas[pheno].url
        file = gwas[pheno].file
        isdir(path) || mkdir(path)
        isfile("$(path)/$(file)") || download(url, "$(path)/$(file)")
    end
end