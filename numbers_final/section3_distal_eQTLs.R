#We asked whether the effects of eQTLs varied across the different stages of the cell-cycle. Of the 2,945 total local eQTLs detected in the three crosses, only 116 (4%) showed 
#significant interactions between the eQTL effect and the cell-cycle stage at an FDR of 5%.

combined_objects$cis_table %>% summarise(has_int=sum(FDR < 0.05),n=n())
combined_objects$cis_table %>% summarise(has_int=sum(has_interaction),n=n())
116/2945

#  In contrast, 790 (24.4%) of 3,238 distal eQTLs showed significant
#interactions with the cell-cycle stage (OR=7.8, Fisher’s exact test, P < 2.2e-16),
combined_objects$hotspot_peaks %>%  summarise(has_int=sum(has_interaction_trans),n=n())


mat1 = matrix(c(116,2829,790,2448),byrow = T,ncol=2,nrow=2)
mat1_res= fisher.test(mat1)
1/mat1_res$estimate
