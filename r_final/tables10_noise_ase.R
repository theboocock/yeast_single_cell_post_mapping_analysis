#out_noise_tmp = combined_objects$noise$ASE %>% inner_join(genes_name_trans,by=c("gene"="gene_id"))

out_noise = pred_df5 %>% dplyr::select(gene.x,gene_name,estimate.cond,p.value.cond,p_adj_ase,theta,p.value.disp,p_adj_disp,good_disp)  %>% mutate(good_disp=!good_disp)#
#%>% ggplot(aes(y=mean_shift,x=estimate.cond)) + geom_point()

#out_noise = out_noise_tmp %>% dplyr::select(gene,gene_name,estimate.cond,p.value.cond,p_adj_ase,estimate.disp,p.value.disp,p_adj_disp)
col_s_out = c("transcript","gene name","estimate (average expression)","p-value (average expression)","adjusted p-value (average expression)",
              "estimate (noise)","p-value (noise)","adjusted p-value (noise)","Overlaps global trend line")
colnames(out_noise)  = col_s_out
c_l = split(out_noise,pred_df5$cross.x)
names(c_l) = c("C","A","B")
c_l = c_l[c("A","B","C")]
#col_s_out 

openxlsx::write.xlsx(c_l,file="tables/s10.xlsx")

