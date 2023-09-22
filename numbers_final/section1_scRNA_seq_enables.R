#We captured a median of 1,514 unique RNA molecules (unique molecular identifiers; henceforth UMIs) 
#and a median of 1,091 expressed SNPs per cell (Table S4).
summary_table[1,]
#Using this classification approach, we found that expression of 2,139 genes displayed
#significant variation by cell cycle stage (likelihood ratio test, false-discovery rate FDR < 0.05; Table S5). 
sum(ap_icc$CC.q < 0.05,na.rm=T)

#We observed a median of 17 cells per segregant, with 277 of the segregants 
#sampled more than ten times (Figure S4)

median(seg_match_df$Freq)
sum(seg_match_df$Freq > 10)
#The genotypes measured from scRNA-seq data were in high agreement with those obtained from whole-genome sequencing 
#of the same strains (median genotype agreement 92.5%). 
median(ap_df$barcode.features$accuracy)
#
#We mapped 770 local eQTLs at an FDR of 5% with the HMM-based genotypes, 
#and 697 with the matched genotypes obtained from whole-genome sequencing 
#of the segregants; 611 eQTLs were detected in both analyses (Table S6). 
nrow(ap_combined %>% filter(p_adj.old < 0.05))
nrow(ap_combined %>% filter(p_adj.hmm < 0.05))
sum(ap_combined$sig_both)
#We further compared the local eQTL effects for all 4,901 tested transcripts, 
#regardless of statistical significance,
#and found that they were highly correlated between the two sets of genotypes
cor.test(ap_combined$beta.hmm,ap_combined$beta.old,method="spearman")
