genofile = snpgdsOpen("data/joseph//1012.gds")



snp_id = read.gdsn(index.gdsn(genofile,"snp.id"))
position = read.gdsn(index.gdsn(genofile,"snp.position"))
chromosome = read.gdsn(index.gdsn(genofile,"snp.chromosome"))



all_snps_gr = data.frame(seqnames=chromosome,start=position,width=1) %>% as_granges() 
gr_hsp12 = as_granges(data.frame(seqnames="VI",start=106323,end=108084))

aa = findOverlaps(gr_hsp12,all_snps_gr)
maf[subjectHits(aa)][12]

all_snps_gr %>% join_overlap_inner(gr_hsp12)


vcf_annotation_16 %>% filter()



byxrm = vcf_annotation_16 %>% filter((grepl("BY",ref_strains) & grepl("RM",alt_strains)) |  (grepl("BY",ref_strains) & grepl("RM",alt_strains)))
byxrm2 = byxrm %>%  filter((grepl("YJM981",ref_strains) & grepl("CBS2888",alt_strains)) |  (grepl("CBS2888",ref_strains) & grepl("YJM981",alt_strains)))

gr_hsp122 = as_granges(data.frame(seqnames="chrVI",start=106000,end=107724))

join_overlap_inner(byxrm2,gr_hsp122)



snpgdsVCF2GDS("/media/theboocock/scratch/vcf_tmp/chrVI.vcf.gz", "data/joseph/indel.gds", method = "copy.num.of.ref")

genofile2 = snpgdsOpen("data/joseph/indel.gds")
#

strains_cross = joseph_annotation$Standardized.name[grep("CBS2888|YJM981|RM",joseph_annotation$Isolate.name)]

strains_cross2 = joseph_annotation$Isolate.name[grep("CBS2888|YJM981|RM",joseph_annotation$Isolate.name)]

snpgdsSummary("data/joseph/indel.gds")
idx_in = grep("CBS2888|YJM981|RM",joseph_annotation$Isolate.name)

sample.int()

samples = read.gdsn(index.gdsn(genofile2,"sample.id"))

idx = which(samples %in% strains_cross)
gt_indel = read.gdsn(index.gdsn(genofile2,"genotype"))

gt_indel[idx,]


#2.32 =
  
prop_with_plasmid = 86000/127e6

1.5/2.32*86000
