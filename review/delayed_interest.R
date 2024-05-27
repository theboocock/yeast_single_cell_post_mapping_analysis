

strains = unique(c("YJM451","YJM1399","YJM1386","YJM554","YJM1339","YJM1418","YJM1252","YJM1444","YJM1338","YJM451"))

grep_str= paste(strains,sep="",collapse = "|")


het_df[grep(grep_str,het_df$Isolate.name),]

genofile2 = snpgdsOpen("data/rare.gds") # MAF > 0.01 & MAF < 0.04??
pos2 = read.gdsn(index.gdsn(genofile2,"snp.position"))
chrom2 =  read.gdsn(index.gdsn(genofile2,"snp.chromosome"))
allele2 =  read.gdsn(index.gdsn(genofile2,"snp.allele"))

idx_rare = which(pos2 == 113512 & chrom2 == "VIII")
gt2 =  read.gdsn(index.gdsn(genofile2,"genotype"))
gt_by = gt2[,idx_rare]
#i#dx_mosaic2 = which(rownames(dissim$diss) %in% m1$standardized_name2)
#snpgdsVCF2GDS("/media/theboocock/scratch/vcf_tmp/maf_0.05_joseph.vcf.gz","1012.gds",method="biallelic.only")
#snpgdsSummary("1012.gds")
#genofile = snpgdsOpen("1012.gds")
#snpgdsClose(genofile)
#dissim = snpgdsDiss(genofile,autosome.only = F)
#rownames(dissim$diss)= dissim$sample.id
#colnames(dissim$diss) = dissim$sample.id

### 

#rep("Mosaic",joseph_annotation$Clades),]
#m1$standardized_name2
by_var  = data.frame(sample=rownames(dissim$diss),gt=gt_by)# %>%


het_df_w_by = het_df %>% inner_join(by_var,by=c("sample"))

het_df_w_by[grep(grep_str,het_df_w_by$Isolate.name),]
