# This locus is also a distal eQTL hotspot that influenced the expression of 51 genes (Figure 2A) #
nrow(cross_data$B$trans$hotspot_enrichments_and_overlaps$`chrVIII:46887-140660_51`$directional_hotspot_distal)
# binomial test

#  We observed that for 36 of the 50 genes affected by the hotspot and detected in our single-cell validation dataset, the sign of the expression difference was consistent between 
#the eQTL effect and the W82R validation experiment (binomial test, p=0.0026; Table S12). I
sum(bl$Beta.y * bl$avg_log2FC < 0,na.rm=T)
binom.test(x=36,n=50)

#  Importantly, the gene expression difference in the W82R experiment was statistically significant and concordant with the eQTL effect for all
#6 mating-related genes affected by this hotspot (AGA1, AGA2, MFA1, STE2, FUS3, PRM5), 
mating_rel = bl[bl$gene_name %in% c("AGA1","AGA2","MFA1","STE2","FUS3","PRM5"),]
mating_rel %>% filter(avg_log2FC * Beta.y < 0) %>% filter(p_val_adj < 0.05)
# Frequency in the different cell-cycle stages ###
#Consistent with the QTL effect, cells with the 82R allele were overrepresented in
#G1 (46.2% in G1 vs. 42.8% in other stages, logistic regression, P<1e-52; Figure S8).
# A  =420/421
# B = 416/417
aaa = table(gpa_c$cell_cycle == "G1",gpa_c$dataset)
total = apply(aaa,1,sum)
aaa/total
out_df_gpa = data.frame()
out_df_gpa1_without_cc = data.frame()
for(cc in unique(gpa_c$cell_cycle)){
  y = as.numeric(gpa_c$cell_cycle == cc)#~ gpa_c$dataset)
  aa = glm(y ~ log(gpa_c$nCount_RNA) + gpa_c$dataset,family="binomial")
  af = (drop1(aa,test="Chisq"))
  aa1 = glm(y ~ gpa_c$dataset,family="binomial")
  aa = glm(y ~ log(gpa_c$nCount_RNA) * gpa_c$dataset)
  aa = glm(y ~  gpa_c$dataset)
  aa = emmeans::emmeans(aa,spec=~dataset)
  p1 = pairs(aa)
  out_df_gpa = rbind(out_df_gpa,data.frame(p1) %>% mutate(lrt=af$`Pr(>Chi)`[3],cell_cycle=cc))
  af = (drop1(aa1,test="Chisq"))
  aa1 = emmeans::emmeans(aa1, spec=~dataset)
  p1=pairs(aa1)
  out_df_gpa1_without_cc = rbind(out_df_gpa1_without_cc,data.frame(p1) %>% mutate(lrt=af$`Pr(>Chi)`[2],cell_cycle=cc))
}

out_df_gpa
out_df_gpa1_without_cc

aa = out_df_gpa[3,]

#We measured growth rates of the engineered strains and found that
#cells with the 82R allele grew slower than those with the 82W allele 
#(relative fitness=0.993,T=-2.592, p=0.0268), but that this effect was smaller
#than that observed for the S469I variant 
#(relative fitness of the 469I allele=0.991,T=-5.291, p<0.001; Figure 5A). 


m1 = (lm((growth_gpa_strains$doubling) ~
           growth_gpa_strains$ID + growth_gpa_strains$plate))

library(emmeans)
em = emmeans::emmeans(m1,specs=~ID)
em_df = em %>% as_data_frame() 

em_df$emmean/em_df$emmean[3]
pairs(em)
.440/.444
#summary(m1)
# The 469I allele is known to improve the efficiency of mating, a difference we successfully replicated (relative 469I mating efficiency=110%, T=9.73, p<0.001). We observed that the 82R allele also increased mating efficiency (relative mating efficiency=108%, T=6.75, p<0.001), but to a lesser extent 
#than the 469I allele (82R mating efficiency compared to 469I=97.2%, T=-2.98, p<0.0001; Figures 5B and S9). 

mm = lm(mating_efficiency_gpa1_strains$norm ~ mating_efficiency_gpa1_strains$strain_new)
aaaa = emmeans(mm,spec=~strain_new)
aaaa
pairs(aaaa)

#
#
# he 82R allele is common (20.5%) (Figure 5
sum(het_df$gpa)/(nrow(het_df) * 2)

het_df %>% group_by(Clades_trim) %>% summarise(gpa_freq=sum(gpa)/(2*n())) %>% arrange(-gpa_freq)

# Peter et al. identified four groups of mosaic strains, which are characterized by admixture of two or more different lineages through outbreeding, and we observed that the 82R allele 
#is enriched in these mosaic strains (allele frequency = 45.3%, permutation test p=0.007).


idx_all = (maf < maf[idx] + .02 & maf > maf[idx] -0.02)
idx_all[idx] =F
idx_all = which(idx_all)

tmp_gt2 = gt3[,idx]
len =  sum(!is.na(tmp_gt2[idx_mosaic]))
mosaic_prop= sum(tmp_gt2[idx_mosaic],na.rm=T)/(len * 2 )
mosaic_prop

props = c()
set.seed(1)
for(idx2 in sample(idx_all,size=1e4,replace = T)){
  tmp_gt2 = gt3[,idx2]
  #i  
  len =  sum(!is.na(tmp_gt2[idx_mosaic]))
  mosaic_prop2 = sum(tmp_gt2[idx_mosaic],na.rm=T)/(len * 2 )
  #mosaic_prop2= ifelse(mosaic_prop2 > 0.5,1-mosaic_prop2,mosaic_prop2)
  props = c(props, mosaic_prop2)
}
sum(mosaic_prop <= props)/10000

# 0.0068


# We compared homozygous (<5% heterozygous sites) and heterozygous (>5% heterozygous sites) strains, as defined by Peter et al., and found that the 
#82R allele is enriched in heterozygous strains (OR=3.3, Fisher’s exact test, p<1e-07) 
#and associated with higher rates of heterozygosity (Wilcoxon rank sum test, p<1e-15; Figure S11). W
ccc = het_df %>% filter(gpa != 1) %>% filter(Ploidy != 1)
fisher.test(table(ccc$gpa,ccc$Zygosity))
1/.300996
summary(lm(ccc$het_sites_prop ~ ccc$gpa))

wilcox.test(ccc$het_sites_prop ~ ccc$gpa)

### Allele-frequencies in the paper
## /media/theboocock/scratch/vcf_tmp ##
#chrVIII	71885
# PHO84 P259L AC=2016;AF=0.994;AN=2022#
# 0.3%
# GPA1 S469I  AC=1984;AF=0.979;AN=2022;BaseQRankSum=1.72;
# 0.187%