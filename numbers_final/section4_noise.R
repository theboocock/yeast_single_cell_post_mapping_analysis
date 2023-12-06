#combined_objects$noise$ASE  = combined_objects$noise$ASE %>% mutate(sig_new_filt =  (  combined_objects$noise$ASE$p_adj_ase > 0.05 | combined_objects$noise$ASE$estimate.cond*combined_objects$noise$ASE$estimate.disp < 0) & p_adj_disp < 0.05)# & p_adj_disp < 1e-3)# %>%
#We obtained a total of 13,973 single-cell transcriptomes from F1 diploids used to generate the segregants for 
#the three crosses (5,890 for cross A, 2,864 for cross B, and 5,219 for cross C; Table S4).
summary_table[c(6,7,10),]
# We found 3,406 genes with allele-specific effects on average expression levels (668 for cross A, 996 for cross B, and 1,742 for cross C; Table S9). These allele-specific effects were well-correlated
#with local eQTL effects from the eQTL mapping experiments described above (Figure S7). 
sum(combined_objects$cis_eqtl_ase_only_with_onepot$p_adj < 0.05,na.rm=T)
combined_objects$cis_eqtl_ase_only_with_onepot %>% group_by(cross) %>% 
  summarise(sig=sum(p_adj < 0.05),n=n())
#We observed 160 genes with significant 
#interactions between allele-specific expression and cell-cycle stage (Table S9).
sum(combined_objects$cis_eqtl_ase_only_with_onepot$has_interaction.x)
# 874
sum(out_noise$`adjusted p-value (dispersion)` < 0.05)
#377
sum(!out_noise$`Overlaps global trend line`)
#cor.test(combined_objects$noise$ASE$estimate.cond,combined_objects$noise$ASE$theta,method="spearman",use = "pairwise.complete.obs")
cor.test(em_subset$theta,em_subset$emmean.mean,method="spearman")
#cor(ab$estimate.cond,ab$estimate.disp,method="spearman",use="pairwise.complete.obs")
#ab$is_si
ab %>% filter(gene == "YFL014W") %>% filter(cross == "A" | cross == "3004")
