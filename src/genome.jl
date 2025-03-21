const GRCh38_highld = [
    (chr = "1", range1 = 47761740, range2 = 51761740),
    (chr = "1", range1 = 125169943, range2 = 125170022),
    (chr = "1", range1 = 144106678, range2 = 144106709),
    (chr = "1", range1 = 181955019, range2 = 181955047),
    (chr = "2", range1 = 85919365, range2 = 100517106),
    (chr = "2", range1 = 135275091, range2 = 135275210),
    (chr = "2", range1 = 182427027, range2 = 189427029),
    (chr = "3", range1 = 47483505, range2 = 49987563),
    (chr = "3", range1 = 83368158, range2 = 86868160),
    (chr = "5", range1 = 44464140, range2 = 51168409),
    (chr = "5", range1 = 129636407, range2 = 132636409),
    (chr = "6", range1 = 25391792, range2 = 33424245),
    (chr = "6", range1 = 57788603, range2 = 58453888),
    (chr = "6", range1 = 61109122, range2 = 61424451),
    (chr = "6", range1 = 139637169, range2 = 142137170),
    (chr = "7", range1 = 54964812, range2 = 66897578),
    (chr = "8", range1 = 8105067, range2 = 12105082),
    (chr = "8", range1 = 43025699, range2 = 48924888),
    (chr = "8", range1 = 110918594, range2 = 113918595),
    (chr = "9", range1 = 64198500, range2 = 64200392),
    (chr = "10", range1 = 36671065, range2 = 43184546),
    (chr = "11", range1 = 88127183, range2 = 91127184),
    (chr = "12", range1 = 32955798, range2 = 41319931),
    (chr = "20", range1 = 33948532, range2 = 36438183),
]

const GRCh37_highld = [
    (chr = "1", range1 = 48000000, range2 = 52000000),
    (chr = "2", range1 = 86000000, range2 = 100500000),
    (chr = "2", range1 = 134500000, range2 = 138000000),
    (chr = "2", range1 = 183000000, range2 = 190000000),
    (chr = "3", range1 = 47500000, range2 = 50000000),
    (chr = "3", range1 = 83500000, range2 = 87000000),
    (chr = "3", range1 = 89000000, range2 = 97500000),
    (chr = "5", range1 = 44500000, range2 = 50500000),
    (chr = "5", range1 = 98000000, range2 = 100500000),
    (chr = "5", range1 = 129000000, range2 = 132000000),
    (chr = "5", range1 = 135500000, range2 = 138500000),
    (chr = "6", range1 = 24000000, range2 = 36000000),
    (chr = "6", range1 = 57000000, range2 = 64000000),
    (chr = "6", range1 = 140000000, range2 = 142500000),
    (chr = "7", range1 = 55000000, range2 = 66000000),
    (chr = "8", range1 = 7000000, range2 = 13000000),
    (chr = "8", range1 = 43000000, range2 = 50000000),
    (chr = "8", range1 = 112000000, range2 = 115000000),
    (chr = "10", range1 = 37000000, range2 = 43000000),
    (chr = "11", range1 = 46000000, range2 = 57000000),
    (chr = "11", range1 = 87500000, range2 = 90500000),
    (chr = "12", range1 = 33000000, range2 = 40000000),
    (chr = "12", range1 = 109500000, range2 = 112000000),
    (chr = "20", range1 = 32000000, range2 = 34500000),
    (chr = "17", range1 = 42900000, range2 = 45100000)
]
# https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)

const GRCh38_totlength = Dict(
    "1" => 248956422,
    "2" => 242193529,
    "3" => 198295559,
    "4" => 190214555,
    "5" => 181538259,
    "6" => 170805979,
    "7" => 159345973,
    "8" => 145138636,
    "9" => 138394717,
    "10" => 133797422,
    "11" => 135086622,
    "12" => 133275309,
    "13" => 114364328,
    "14" => 107043718,
    "15" => 101991189,
    "16" => 90338345,
    "17" => 83257441,
    "18" => 80373285,
    "19" => 58617616,
    "20" => 64444167,
    "21" => 46709983,
    "22" => 50818468,
    "X" => 156040895,
    "Y" => 57227415,
)

const GRCh37_totlength = Dict(
    "1" => 249250621,
    "2" => 243199373,
    "3" => 198022430,
    "4" => 191154276,
    "5" => 180915260,
    "6" => 171115067,
    "7" => 159138663,
    "8" => 146364022,
    "9" => 141213431,
    "10" => 135534747,
    "11" => 135006516,
    "12" => 133851895,
    "13" => 115169878,
    "14" => 107349540,
    "15" => 102531392,
    "16" => 90354753,
    "17" => 81195210,
    "18" => 78077248,
    "19" => 59128983,
    "20" => 63025520,
    "21" => 48129895,
    "22" => 51304566,
    "X" => 155270560,
    "Y" => 59373566,
)
# https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
