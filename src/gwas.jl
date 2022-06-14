gwas = Dict(
    # Psychiatric
    "asd" => (url = "https://figshare.com/ndownloader/files/28169292",
        pmid = "30804558",
        title = "Autism (Grove et al. 2019)",
        file = "iPSYCH-PGC_ASD_Nov2017.gz"),
    "scz" => (url = "https://figshare.com/ndownloader/files/28169757",
        pmid = "35396580",
        title = "Schizophrenia (PGC3)",
        file = "PGC3_SCZ_wave3_public.v2.tsv.gz"),
    "bd" => (url = "https://figshare.com/ndownloader/files/26603681",
        pmid = "34002096",
        title = "Bipolar (Mullins et al. 2021)",
        file = "pgc-bip2021-all.vcf.tsv.gz"),
    "sczvsbd" => (url = "https://figshare.com/ndownloader/files/28169358",
        pmid = "29906448",
        title = "SCZ vs BD (Ruderfer et al. 2018)",
        file = "SCZvsBD.sumstats.gz"),
    "adhd" => (url = "https://figshare.com/ndownloader/files/28169253",
        pmid = "30478444",
        title = "ADHD (Demontis et al. 2019)",
        file = "daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz"),
    "anorexia" => (url = "https://figshare.com/ndownloader/files/28169271",
        pmid = "31308545",
        title = "Anorexia (Watson et al. 2019)",
        file = "pgcAN2.2019-07.vcf.tsv.gz"),
    "cross" => (url = "https://figshare.com/ndownloader/files/28169382",
        pmid = "31835028",
        title = "Cross disorder (PGC 2019)",
        file = "pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt.gz"),
    "mdd" => (url = "https://datashare.ed.ac.uk/bitstream/handle/10283/3203/PGC_UKB_depression_genome-wide.txt",
        pmid = "30718901",
        title = "Major depression (Howard et al. 2019)",
        file = "PGC_UKB_depression_genome-wide.txt"),
    # Neurologic
    "alz" => (url = "https://ctg.cncr.nl/documents/p1651/AD_sumstats_Jansenetal_2019sept.txt.gz",
        pmid = "30617256",
        title = "Alzheimer disease (Jansen et al. 2019)",
        file = "AD_sumstats_Jansenetal_2019sept.txt.gz"),
    "als" => (url = "https://surfdrive.surf.nl/files/index.php/s/E5RetKw10hC3jXy/download?path=%2F&files=ALS_GWAS_RVB_SMR_sumstats.tar.gz",
        pmid = "34873335",
        title = "ALS (van Rheenen et al. 2021)",
        file = "ALS_sumstats_EUR_only.txt.gz"),
    "parkinson" => (url = "https://drive.google.com/file/d/1FZ9UL99LAqyWnyNBxxlx6qOUlfAnublN",
        pmid = "31701892",
        title = "Parkinson's disease (Nalls et al. 2019)",
        file = "nallsEtAl2019_excluding23andMe_allVariants.tab"),
    "stroke" => (url = "https://www.kp4cd.org/datasets/stroke",
        pmid = "29531354",
        title = "Stroke (Malik et al. 2018)",
        file = "metastroke.all.chr.bp"),
    # Neurocognitive
    "left-hand" => (url = "https://evansgroup.di.uq.edu.au/GWAS_RESULTS/HANDEDNESS/LeftHandedness_MetaAnalysis_UKBB_IHC.txt.gz",
        pmid = "32989287",
        title = "Left-handedness (Cuellar-Partida et al. 2020)",
        file = "LeftHandedness_MetaAnalysis_UKBB_IHC.txt.gz"),
    "cp" => (url = "https://thessgac.com/papers/",
        pmid = "30038396",
        title = "Cognitive performance (Lee et al. 2018)",
        file = "GWAS_CP_all.txt"),
    "ea" => (url = "https://thessgac.com/papers/",
        pmid = "35361970",
        title = "Educational attainment (Okbay et al. 2022)",
        file = "EA4_additive_excl_23andMe.txt.gz"),
    "birth" => (url = "https://thessgac.com/papers/",
        pmid = "27798627",
        title = "Age at first birth (Barban et al. 2016)",
        file = "AgeFirstBirth_Pooled.txt"),
    "children" => (url = "https://thessgac.com/papers/",
        pmid = "27798627",
        title = "Number of children (Barban et al. 2016)",
        file = "NumberChildrenEverBorn_Pooled.txt"),
    "risk" => (url = "https://thessgac.com/papers/",
        pmid = "30643258",
        title = "Risk tolerance (Linner et al. 2019)",
        file = "RISK_GWAS_MA_UKB+replication.txt"),
    "neuroticism" => (url = "https://ctg.cncr.nl/documents/p1651/sumstats_neuroticism_ctg_format.txt.gz",
        pmid = "29942085",
        title = "Neuroticism (Nagel et al. 2018)",
        file = "sumstats_neuroticism_ctg_format.txt.gz"),
    "smoking" => (url = "https://conservancy.umn.edu/bitstream/handle/11299/201564/CigarettesPerDay.txt.gz",
        pmid = "30643251",
        title = "Cigarettes per day (Liu et al. 2019)",
        file = "CigarettesPerDay.txt.gz"),
    "alcohol" => (url = "https://conservancy.umn.edu/bitstream/handle/11299/201564/DrinksPerWeek.txt.gz",
        pmid = "30643251",
        title = "Drinks per week (Liu et al. 2019)",
        file = "DrinksPerWeek.txt.gz"),
    "intelligence" => (url = "https://ctg.cncr.nl/documents/p1651/SavageJansen_IntMeta_sumstats.zip",
        pmid = "29942086",
        title = "Intelligence (Savage et al. 2018)",
        file = "SavageJansen_2018_intelligence_metaanalysis.txt"),
    "insomnia" => (url = "https://ctg.cncr.nl/documents/p1651/Insomnia_sumstats_Jansenetal.txt.gz",
        pmid = "30804565",
        title = "Insomnia (Jansen et al. 2019)",
        file = "Insomnia_sumstats_Jansenetal.txt.gz"),
    "icv" => (url = "https://ctg.cncr.nl/documents/p1651/meta_analysis_BV_Jansenetal_2020.sumstats.txt.gz",
        pmid = "33154357",
        title = "Intracranial volume (Jansen et al. 2020)",
        file = "meta_analysis_BV_Jansenetal_2020.sumstats.txt.gz"),
    # Ophtho
    "amd" => (url = "http://csg.sph.umich.edu/abecasis/public/amd2015/Fritsche_2015_AdvancedAMD.txt.gz",
        pmid = "26691988",
        title = "AMD (Fritsche et al. 2016)",
        file = "Fritsche_2015_AdvancedAMD.txt.gz"),
    "cataract" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90014001-GCST90015000/GCST90014268/GCST90014268_buildGRCh37.tsv",
        pmid = "34127677",
        title = "Cataract (Choquet et al. 2021)",
        file = "GCST90014268_buildGRCh37.tsv"),
    "glaucoma" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009722/MTAG_glaucoma_four_traits_summary_statistics.txt",
        pmid = "31959993",
        title = "Galucoma (Craig et al. 2020)",
        file = "MTAG_glaucoma_four_traits_summary_statistics.txt"),
    # Gyn
    "menarche" => (url = "https://www.reprogen.org/Menarche_1KG_NatGen2017_WebsiteUpload.zip",
        pmid = "28436984",
        title = "Age at menarche (Day et al. 2017)",
        file = "Menarche_1KG_NatGen2017_WebsiteUpload.txt"),
    "menopause" => (url = "https://www.reprogen.org/reprogen_ANM_201K_170621.txt.gz",
        pmid = "34349265",
        title = "Age at menoapuse (Ruth et al. 2021)",
        file = "reprogen_ANM_201K_170621.txt.gz"),
    "pcos" => (url = "https://www.repository.cam.ac.uk/bitstream/handle/1810/283491/PCOS_summary_data_19092018.tar.gz",
        pmid = "30566500",
        title = "PCOS (Day et al. 2018)",
        file = "PCOS_summary_data_19092018.txt"),
    # Autoimmune
    "lupus" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003156/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz",
        pmid = "26502338",
        title = "Lupus (Bentham et al. 2015)",
        file = "bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz"),
    "ra" => (url = "http://plaza.umin.ac.jp/~yokada/datasource/files/GWASMetaResults/RA_GWASmeta_European_v2.txt.gz",
        pmid = "24390342",
        title = "Rheumatoid arthritis (Okada et al. 2014)",
        file = "RA_GWASmeta_European_v2.txt.gz"),
    "t1d" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000529/harmonised/33830302-GCST90000529-EFO_0001359-Build37.f.tsv.gz",
        pmid = "33830302",
        title = "Type I diabetes (Inshaw et al. 2021)",
        file = "33830302-GCST90000529-EFO_0001359-Build37.f.tsv.gz"),
    "celiac" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005523/harmonised/22057235-GCST005523-EFO_0001060-Build37.f.tsv.gz",
        pmid = "22057235",
        title = "Celiac disease (Trynka et al. 2011)",
        file = "22057235-GCST005523-EFO_0001060-Build37.f.tsv.gz"),
    "uc" => (url = "ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/human/2016-11-07/",
        pmid = "28067908",
        title = "Ulcerative colitis (de Lange et al. 2017)",
        file = "uc_build37_45975_20161107.txt.gz"),
    "crohns" => (url = "ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/human/2016-11-07/",
        pmid = "28067908",
        title = "Crohn's disease (de Lange et al. 2017)",
        file = "cd_build37_40266_20161107.txt.gz"),
    "ibd" => (url = "ftp://ftp.sanger.ac.uk/pub/project/humgen/summary_statistics/human/2016-11-07/",
        pmid = "28067908",
        title = "IBD (de Lange et al. 2017)",
        file = "ibd_build37_59957_20161107.txt.gz"),
    "as" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005529/harmonised/23749187-GCST005529-EFO_0003898-Build37.f.tsv.gz",
        pmid = "23749187", 
        title = "Ankylosing spondylitis (Cortes et al. 2013)",
        file = "23749187-GCST005529-EFO_0003898-Build37.f.tsv.gz"),
    "pbc" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003129/harmonised/26394269-GCST003129-EFO_1001486-Build37.f.tsv.gz",
        pmid = "26394269",
        title = "Primary biliary cirrhosis (Cordell et al. 2015)",
        file = "26394269-GCST003129-EFO_1001486-Build37.f.tsv.gz"),
    # Respiratory
    "asthma" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010042/HanY_prePMID_asthma_UKBB.txt.gz",
        pmid = "32296059",
        title = "Asthma (Han et al. 2020)",
        file = "HanY_prePMID_asthma_UKBB.txt.gz"),
    # MSK
    "oa" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005811/Zengini_29559693_SR.txt.gz",
        pmid = "29559693",
        title = "Osteoarthritis (Zengini et al. 2018)",
        file = "Zengini_29559693_SR.txt.gz"),
    # Digestive
    "gerd" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000514/GCST90000514_buildGRCh37.tsv.gz",
        pmid = "34187846",
        title = "GERD (Ong et al. 2022)",
        file = "GCST90000514_buildGRCh37.tsv.gz"),
    "barrett" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90000001-GCST90001000/GCST90000515/GCST90000515_buildGRCh37.tsv.gz", 
        pmid = "34187846",
        title = "Barrett's esophagus (Ong et al. 2022)",
        file = "GCST90000515_buildGRCh37.tsv.gz"),
    # Anthropometric
    "height" => (url = "https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz",
        pmid = "30124842",
        title = "Height (Yengo et al. 2018)",
        file = "Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz"),
    "weight" => (url = "https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz",
        pmid = "30124842",
        title = "Weight (Yengo et al. 2018)",
        file = "Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz"),
    # Hormone
    "tsh" => (url = "http://csg.sph.umich.edu/willer/public/TSH2020/2020_Zhou_et_al_TSH_meta_HUNT-MGI-ThyroidOmics.txt",
        pmid = "32769997",
        title = "TSH level (Zhou et al. 2020)",
        file = "2020_Zhou_et_al_TSH_meta_HUNT-MGI-ThyroidOmics.txt"),
    # Cardiovascular
    "cad" => (url = "http://www.cardiogramplusc4d.org/data-downloads/",
        pmid = "26343387",
        title = "CAD (Nikpay et al. 2015)",
        file = "cad.add.160614.website.txt"),
    "afib" => (url = "https://personal.broadinstitute.org/ryank/AF_HRC_GWAS_ALLv11.zip",
        pmid = "29892015",
        title = "Atrial fibrillation (Roselli et al. 2018)",
        file = "AF_HRC_GWAS_ALLv11.txt"),
    "dbp" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006630/Evangelou_30224653_DBP.txt.gz", 
        pmid = "30224653",
        title = "DBP (Evangelou et al. 2018)",
        file = "Evangelou_30224653_DBP.txt.gz"),
    "sbp" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz",
        pmid = "30224653",
        title = "SBP (Evangelou et al. 2018)",
        file = "Evangelou_30224653_SBP.txt.gz"),
    # Kidney
    "ckd" => (url = "http://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
        pmid = "31152163",
        title = "CKD (Wuttke et al. 2019)",
        file = "CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz"),
    # Endocrine
    "t2d" => (url = "https://diagram-consortium.org/downloads.html",
        pmid = "30297969",
        title = "Type II diabetes (Mahajan et al. 2018)",
        file = "Mahajan.NatGenet2018b.T2D.European.txt"),
    # Respiratory
    "covid" => (url = "https://storage.googleapis.com/covid19-hg-public/20210415/results/20210607/COVID19_HGI_A2_ALL_leave_23andme_20210607.b37.txt.gz",
        pmid = "34237774",
        title = "COVID-19 (2021)",
        file = "COVID19_HGI_A2_ALL_leave_23andme_20210607.b37.txt.gz"),
    # Derm
    "acne" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90091001-GCST90092000/GCST90092000/GCST90092000_buildGRCh37.tsv",
        pmid = "35132056",
        title = "Acne (Mitchell et al. 2022)",
        file = "GCST90092000_buildGRCh37.tsv.gz"),
    "bald" => (url = "http://www.psy.ed.ac.uk/ccace/downloads/Hagenaars2017_UKB_MPB_summary_results.zip",
        pmid = "28196072",
        title = "Male-pattern baldness (Hagenaars et al. 2017)",
        file = "Hagenaars2017_UKB_MPB_summary_results_10Feb2016.txt"),
    "ad" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003184/EAGLE_AD_no23andme_results_29072015.txt",
        pmid = "26482879",
        title = "Atopic dermatitis (Paternoster et al. 2015)",
        file = "EAGLE_AD_no23andme_results_29072015.txt"),
    "vitiligo" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004785/",
        pmid = "27723757",
        title = "Vitiligo (Jin et al. 2016)",
        file = "GWAS123.txt.gz"),
    "psoriasis" => (url = "https://gwas.mrcieu.ac.uk/files/ukb-b-10537/ukb-b-10537.vcf.gz",
        pmid = "",
        title = "Psoriasis (UK Biobank)",
        file = "ukb-b-10537.vcf.gz"),
    "bcc" => (url = "https://gwas.mrcieu.ac.uk/files/ukb-b-8837/ukb-b-8837.vcf.gz",
        pmid = "",
        title = "Basal cell carcinoma (UK Biobank)",
        file = "ukb-b-8837.vcf.gz"),
    "melanoma" => (url = "https://gwas.mrcieu.ac.uk/files/ukb-b-12915/ukb-b-12915.vcf.gz",
        pmid = "",
        title = "Melanoma (UK Biobank)",
        file = "ukb-b-12915.vcf.gz"),
    # Cancer
    "breast" => (url = "https://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/oncoarray-and-combined-summary-result/gwas-summary-associations-breast-cancer-risk-2020/",
        pmid = "32424353",
        title = "Breast cancer (Zhang et al. 2020)",
        file = "icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt"),
    "prostate" => (url = "https://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/oncoarray-and-combined-summary-result/gwas-summary-associations-breast-cancer-risk-2020/",
        pmid = "23065704",
        title = "Prostate cancer (Schumacher et al. 2018)",
        file = "meta_v3_onco_euro_overall_ChrAll_1_release.txt"),
    "ec" => (url = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006464/GCST006464_GRCh37.tsv.gz",
        pmid = "30093612",
        title = "Endometrial cancer (O'Mara et al. 2018)",
        file = "GCST006464_GRCh37.tsv.gz"),
    # CBC
    "rbc" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_RBC_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Red blood cell count (Chen et al. 2020)",
        file = "BCX2_RBC_EA_GWAMA.out.gz"),
    "hgb" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_HGB_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Hemoglobin concentration (Chen et al. 2020)",
        file = "BCX2_HGB_EA_GWAMA.out.gz"),
    "hct" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_HCT_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Hematocrit (Chen et al. 2020)",
        file = "BCX2_HCT_EA_GWAMA.out.gz"),
    "mch" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_MCH_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Mean corpuscular hemoglobin (Chen et al. 2020)",
        file = "BCX2_MCH_EA_GWAMA.out.gz"),
    "mcv" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_MCV_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Mean corpuscular volume (Chen et al. 2020)",
        file = "BCX2_MCV_EA_GWAMA.out.gz"),
    "mchc" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_MCHC_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Mean corpuscular hemoglobin concentration (Chen et al. 2020)",
        file = "BCX2_MCHC_EA_GWAMA.out.gz"),
    "rdw" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_RDW_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "RBC distribution width (Chen et al. 2020)",
        file = "BCX2_RDW_EA_GWAMA.out.gz"),
    "wbc" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_WBC_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "White blood cell count (Chen et al. 2020)",
        file = "BCX2_WBC_EA_GWAMA.out.gz"),
    "neu" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_NEU_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Neutrophil count (Chen et al. 2020)",
        file = "BCX2_NEU_EA_GWAMA.out.gz"),
    "lymph" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_LYM_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Lymphocyte count (Chen et al. 2020)",
        file = "BCX2_LYM_EA_GWAMA.out.gz"),
    "mono" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_MON_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Monocyte count (Chen et al. 2020)",
        file = "BCX2_MON_EA_GWAMA.out.gz"),
    "baso" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_BAS_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Basophil count (Chen et al. 2020)",
        file = "BCX2_BAS_EA_GWAMA.out.gz"),
    "eos" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_EOS_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Eosinophil count (Chen et al. 2020)",
        file = "BCX2_EOS_EA_GWAMA.out.gz"),
    "plt" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_PLT_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Platelet count (Chen et al. 2020)",
        file = "BCX2_PLT_EA_GWAMA.out.gz"),
    "mpv" => (url = "http://www.mhi-humangenetics.org/dataset/BCX2_MPV_EA_GWAMA.out.gz",
        pmid = "32888493",
        title = "Mean platelet volume (Chen et al. 2020)",
        file = "BCX2_MPV_EA_GWAMA.out.gz"),
    # Serum/urine biomarker
    "alt" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "ALT (Sinnott-Armstrong et al. 2021)",
        file = "Alanine_aminotransferase.imp.gz"),
    "ast" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "AST (Sinnott-Armstrong et al. 2021)",
        file = "Aspartate_aminotransferase.imp.gz"),
    "astalt" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "AST / ALT (Sinnott-Armstrong et al. 2021)",
        file = "AST_ALT_ratio.imp.gz"),
    "ap" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Alkaline phosphatase (Sinnott-Armstrong et al. 2021)",
        file = "Alkaline_phosphatase.imp.gz"),
    "apoa" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "APOA (Sinnott-Armstrong et al. 2021)",
        file = "Apolipoprotein_A.imp.gz"),
    "apob" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "APOB (Sinnott-Armstrong et al. 2021)",
        file = "Apolipoprotein_B_adjstatins.imp.gz"),
    "crp" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "CRP (Sinnott-Armstrong et al. 2021)",
        file = "C_reactive_protein.imp.gz"),
    "calcium" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Ca²⁺ (Sinnott-Armstrong et al. 2021)",
        file = "Calcium.imp.gz"),
    "cholesterol" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Cholesterol (Sinnott-Armstrong et al. 2021)",
        file = "Cholesterol_adjstatins.imp.gz"),
    "creatinine" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Creatinine (Sinnott-Armstrong et al. 2021)",
        file = "Creatinine.imp.gz"),
    "ucreatinine" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Urine creatinine (Sinnott-Armstrong et al. 2021)",
        file = "Creatinine_in_urine.imp.gz"),
    "cystatinc" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Cystatin C (Sinnott-Armstrong et al. 2021)",
        file = "Cystatin_C.imp.gz"),
    "albumin" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Albumin (Sinnott-Armstrong et al. 2021)",
        file = "Albumin.imp.gz"),
    "dirbilirubin" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Direct bilirubin (Sinnott-Armstrong et al. 2021)",
        file = "Direct_bilirubin.imp.gz"),
    "egfr" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "eGFR (Sinnott-Armstrong et al. 2021)",
        file = "eGFR.imp.gz"),
    "ggt" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "GGT (Sinnott-Armstrong et al. 2021)",
        file = "Gamma_glutamyltransferase.imp.gz"),
    "glucose" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Glucose (Sinnott-Armstrong et al. 2021)",
        file = "Glucose.imp.gz"),
    "hba1c" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "HbA1c (Sinnott-Armstrong et al. 2021)",
        file = "Glycated_haemoglobin_HbA1c.imp.gz"),
    "hdl" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "HDL cholesterol (Sinnott-Armstrong et al. 2021)",
        file = "HDL_cholesterol.imp.gz"),
    "igf1" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "IGF-1 (Sinnott-Armstrong et al. 2021)",
        file = "IGF_1.imp.gz"),
    "ldl" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "LDL cholesterol (Sinnott-Armstrong et al. 2021)",
        file = "LDL_direct_adjstatins.imp.gz"),
    "lpa" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "LPA (Sinnott-Armstrong et al. 2021)",
        file = "Lipoprotein_A.imp.gz"),
    "ualbumin" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Urine microalbumin (Sinnott-Armstrong et al. 2021)",
        file = "Microalbumin_in_urine.imp.gz"),
    "nap" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Non-albumin protein (Sinnott-Armstrong et al. 2021)",
        file = "Non_albumin_protein.imp.gz"),
    "phosphate" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Phosphate (Sinnott-Armstrong et al. 2021)",
        file = "Phosphate.imp.gz"),
    "upotassium" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Urine K⁺ (Sinnott-Armstrong et al. 2021)",
        file = "Potassium_in_urine.imp.gz"),
    "shbg" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "SHBG (Sinnott-Armstrong et al. 2021)",
        file = "SHBG.imp.gz"),
    "usodium" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Urine Na⁺ (Sinnott-Armstrong et al. 2021)",
        file = "Sodium_in_urine.imp.gz"),
    "testosterone" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Testosterone (Sinnott-Armstrong et al. 2021)",
        file = "Testosterone.imp.gz"),
    "totbilirubin" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Total bilirubin (Sinnott-Armstrong et al. 2021)",
        file = "Total_bilirubin.imp.gz"),
    "totprotein" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Total protein (Sinnott-Armstrong et al. 2021)",
        file = "Total_protein.imp.gz"),
    "triglyceride" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Triglyceride (Sinnott-Armstrong et al. 2021)",
        file = "Triglycerides.imp.gz"),
    "urate" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Urate (Sinnott-Armstrong et al. 2021)",
        file = "Urate.imp.gz"),
    "bun" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "BUN (Sinnott-Armstrong et al. 2021)",
        file = "Urea.imp.gz"),
    "vitamind" => (url = "https://github.com/rivas-lab/biomarkers/tree/master/meta_flipfix/figshare_submission",
        pmid = "33462484",
        title = "Vitamin D (Sinnott-Armstrong et al. 2021)",
        file = "Vitamin_D.imp.gz")
    )