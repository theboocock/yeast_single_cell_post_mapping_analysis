combined_objects$noise$ASE = combined_objects$noise$ASE %>% inner_join(genes_name_trans,by=c("gene"="gene_id"))

out_noise = combined_objects$noise$ASE %>% dplyr::select(gene,gene_name,estimate.cond,p.value.cond,p_adj_ase,estimate.disp,p.value.disp,p_adj_disp)
# flip the dispersion estimate # 
out_noise$estimate.disp = out_noise$estimate.disp

#col_s_out 
col_s_out = c("transcript","gene name","estimate (average expression)","p-value (average expression)","adjusted p-value (average expression)",
              "estimate (dispersion)","p-value (dispersion)","adjusted p-value (dispersion)")
colnames(out_noise)  = col_s_out
openxlsx::write.xlsx(out_noise,file="tables/s10.xlsx")


combined_objects$noise$ASE[combined_objects$noise$ASE$sig_new_filt,] %>% group_by(cross) %>% summarise(n=n()) #%>% filter(gene_name == "HSP12")


combined_objects$noise$AScombined_objects$noise$ASE$estimate.disp



combined_objects$noise$ASE_EMMEANS

