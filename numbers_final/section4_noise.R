#combined_objects$noise$ASE  = combined_objects$noise$ASE %>% mutate(sig_new_filt =  (  combined_objects$noise$ASE$p_adj_ase > 0.05 | combined_objects$noise$ASE$estimate.cond*combined_objects$noise$ASE$estimate.disp < 0) & p_adj_disp < 0.05)# & p_adj_disp < 1e-3)# %>%
#We obtained a total of 13,973 single-cell transcriptomes from F1 diploids used to generate the segregants for 
#the three crosses (5,890 for cross A, 2,864 for cross B, and 5,219 for cross C; Table S4).
summary_table[5:7,]

# We found 3,406 genes with allele-specific effects on average expression levels (668 for cross A, 996 for cross B, and 1,742 for cross C; Table S9). These allele-specific effects were well-correlated
#with local eQTL effects from the eQTL mapping experiments described above (Figure S7). 

sum(combined_objects$cis_eqtl_ase_only_with_onepot$p_adj < 0.05,na.rm=T)
combined_objects$cis_eqtl_ase_only_with_onepot %>% group_by(cross) %>% 
  summarise(sig=sum(p_adj < 0.05),n=n())

#We observed 160 genes with significant 
#interactions between allele-specific expression and cell-cycle stage (Table S9).

sum(combined_objects$cis_eqtl_ase_only_with_onepot$has_interaction.x)
# 874
combined_objects$noise$ASE %>% filter(p_adj_disp < 0.05) %>% filter(sig_new_filt)
# 481