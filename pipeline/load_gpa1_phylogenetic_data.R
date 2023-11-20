library(SNPRelate)
genofile = snpgdsOpen("data/joseph//1012.gds")
#snpgdsClose(genofile)
dissim = snpgdsDiss(genofile,autosome.only = F)
rownames(dissim$diss)= dissim$sample.id
colnames(dissim$diss) = dissim$sample.id
joseph_annotation = xlsx::read.xlsx("data/joseph/full_annotations.xlsx",sheetName="peter2018")
joseph_annotation$Clades_trim = (str_trim(joseph_annotation$Clades))
joseph_annotation$mosaic = grepl("Mosaic",joseph_annotation$Clades)
m1 = joseph_annotation[joseph_annotation$mosaic,]
###
idx_mosaic = which(rownames(dissim$diss) %in% m1$standardized_name2)
tree = ape::bionj(dissim$diss)
mid = phangorn::midpoint(tree)
pos = read.gdsn(index.gdsn(genofile,"snp.position"))
chrom =  read.gdsn(index.gdsn(genofile,"snp.chromosome"))
allele =  read.gdsn(index.gdsn(genofile,"snp.allele"))
#Gpa1
idx = which(chrom == "VIII" & pos == "114674")
gt =  read.gdsn(index.gdsn(genofile,"genotype"))
w82r_gt = (2-(gt[,idx]))
gt[gt==3] = NA
gt3 = 2-gt
ac = colSums(gt3,na.rm=T)
tc = colSums(!is.na(gt3))
af = ac / (tc*2)
gt3[,af > 0.5] = 2- gt3[,af  > .5]
maf = ifelse(af > 0.5,1-af,af)
mosaic_prop = sum(w82r_gt[idx_mosaic])/(length(idx_mosaic) * 2 )

het_sites = rowSums(gt == 1,na.rm=T)
het_sites_divisor = rowSums(!is.na(gt))*2
het_sites_prop = het_sites/het_sites_divisor

het_df = data.frame(het_sites_prop=het_sites_prop,gpa=w82r_gt)
het_df$sample = dissim$sample.id
het_df = het_df %>% inner_join(joseph_annotation,by=c("sample"="standardized_name2"))
het_df_filt = het_df %>% filter(gpa != 1) %>% filter(Ploidy != 1)
fisher.test(table(het_df_filt$gpa,het_df_filt$Zygosity))
wilcox.test(het_df_filt$het_sites_prop ~ het_df_filt$gpa)
