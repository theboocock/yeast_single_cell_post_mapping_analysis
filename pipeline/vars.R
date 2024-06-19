sets=list(
  '3004'=c(9,10),
  'A'=c(3,4),
  'B'=c(1,2,11,12),
  'Ap'=c(5,6)
)

number_set = basename(list.dirs("data/out/",recursive = F))

nums = unlist(lapply(str_split(number_set,"_"), function(x){x[1]}))

as.numeric(nums)

#ames(cList)[c(9,10)]

load("data/rr/seg.recoded.RData")
load("data/rr/extracted_average_phenotypes.RData")

cList=list(
  '01_2444_44_1-2'='B',
  '02_2444_44-'='B',
  '03_ByxRM_51-_1-2'='A',
  '04_ByxRM_51-'='A',
  '00_BYxRM_480MatA_1'='A',
  '00_BYxRM_480MatA_2'='A',
  '08_2444_cross_10k_Feb_21'='B',
  '09_2444_cross_5k_Feb_21'='B',
  '10_3004_cross_10k_Feb_21'='3004',
  '11_3004_cross_5k_Feb_21'='3004',
  '07_2444-cross-1'='B',
  '07_2444-cross-2'='B',
  '19_393_10k_May10'='393',
  '20_393_20k_May10'='393',
  '21_3004_10k_May10'='3004',
  '22_3004_20k_May10'='3004'
)
cc.big.table = readr::read_delim("data//cell_cycle_feb02162022.tsv",delim = "\t")
hList=list(
  '01_2444_44_1-2'=2^6.1,
  '02_2444_44-'= 2^7.5,
  '03_ByxRM_51-_1-2'=2^7,
  '04_ByxRM_51-'=2^7.4,
  '00_BYxRM_480MatA_1'=2^6.5,
  '00_BYxRM_480MatA_2'=2^6.5,
  '08_2444_cross_10k_Feb_21'=2^5.2,
  '09_2444_cross_5k_Feb_21'= 2^4.9,
  '10_3004_cross_10k_Feb_21'= 2^5.6,
  '11_3004_cross_5k_Feb_21'=  2^4.9,
  '07_2444-cross-1'  = 2^6,
  '07_2444-cross-2'  = 2^5.5,
  '19_393_10k_May10' = 2^5.2,
  '20_393_20k_May10' = 2^5.1,
  '21_3004_10k_May10'= 2^6,
  '22_3004_20k_May10'= 2^6.7
)


experiments=names(cList)

crosses.to.parents=list(
  '375'=c("M22", "BYa"),           #1
  'A'  =c("BYa", "RMx"),           #2
  '376'=c("RMx", "YPS163a"),       #3
  'B'  =c("YPS163a", "YJM145x"),   #4
  '377'=c("YJM145x", "CLIB413a"),  #5
  '393'=c("CLIB413a", "YJM978x"),  #6
  '381'=c("YJM978x", "YJM454a"),   #7
  '3008'=c("YJM454a", "YPS1009x"),  #8
  '2999'=c("YPS1009x", "I14a"),     #9
  '3000'=c("I14a", "Y10x"),         #10
  '3001'=c("Y10x", "PW5a"),         #11
  '3049'=c("PW5a", "273614xa"),     #12
  '3003'=c("273614xa", "YJM981x"),  #13
  '3004'=c("YJM981x", "CBS2888a"),  #14
  '3043'=c("CBS2888a", "CLIB219x"), #15
  '3028'=c("CLIB219x", "M22")       #16
)

with_ase_rep  = T
if(!with_ase_rep){
  good.dips=list('13_Group1_diploids_3004_2444_3051_7point5K_Feb_21'=names(crosses.to.parents)[c(4,14)],
  '18_3051_May10'=names(crosses.to.parents)[2])
}else{good.dips=list(
  '13_Group1_diploids_3004_2444_3051_7point5K_Feb_21'=names(crosses.to.parents)[c(2,4,14)],
  '14_GP1_3004_3051_274_375_May10'=names(crosses.to.parents)[c(4,14)],
  '18_3051_May10'=names(crosses.to.parents)[2])
}
chrom_lengths = data.frame(chrom=names(BSgenome.Scerevisiae.UCSC.sacCer3), lengths=lengths(BSgenome.Scerevisiae.UCSC.sacCer3))
chrom_lengths_plot = data.frame(chrom=names(BSgenome.Scerevisiae.UCSC.sacCer3)[1:16], pos=lengths(BSgenome.Scerevisiae.UCSC.sacCer3)[1:16])
chrom_lengths_plot = rbind(chrom_lengths_plot,data.frame(chrom=names(BSgenome.Scerevisiae.UCSC.sacCer3)[1:16], pos=0))
chrom_lengths_plot$tchrom_short_f = convert_chrom_to_simple_factor(chrom_lengths_plot$chrom,reverse = T)
chrom_lengths_plot$chrom_short_f = convert_chrom_to_simple_factor(chrom_lengths_plot$chrom)


source("pipeline/load_cc.R")
