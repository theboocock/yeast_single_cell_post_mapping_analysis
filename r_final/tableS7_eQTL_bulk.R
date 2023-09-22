one_pot_bulk_out = one_pot_bulk

one_pot_bulk_out = one_pot_bulk_out  %>% dplyr::select(transcript,closest.marker,Beta,LOD,FDR,coefficient,p.value,p_adj,has_interaction)

#one_pot_bulk_out$LOD = (10^-(one_pot_bulk_out$LOD))
#one_pot_bulk_out
one_pot_bulk_out=one_pot_bulk_out[!is.na(one_pot_bulk_out$closest.marker),]
one_pot_bulk_out = one_pot_bulk_out %>% inner_join(genes_name_trans,by=c("transcript"="gene_id"))
one_pot_bulk_out = one_pot_bulk_out  %>% dplyr::select(transcript,gene_name,closest.marker,Beta,LOD,FDR,coefficient,p.value,p_adj,has_interaction)
one_pot_bulk_out


colnames(one_pot_bulk_out) = c("transcript","gene name","closest marker","estimate (one-pot)","LOD (one-pot)","FDR (one-pot)","estimate (bulk)","p-value (bulk)","adjusted p-value (bulk)","has cell-cycle interaction")
#one_pot_bulk_out %>% write_tsv(file="tables/s7.tsv")
### Get the other two crosses ###


B = cross_data$B$cis$cis_test_with_disp_expr
C = cross_data$`3004`$cis$cis_test_with_disp_expr

B = B  %>% inner_join(genes_name_trans,by=c("transcript"="gene_id"))
C = C  %>% inner_join(genes_name_trans,by=c("transcript"="gene_id"))
# TODO: diploids #
#B$LOD = 10^(-B$LOD)
#C$LOD = 10^(-C$LOD)
col_names_out = c("transcript","gene name","closest marker","estimate (one-pot)","LOD (one-pot)","FDR (one-pot)","has cell-cycle interaction")
B = B %>% dplyr::select(transcript,gene_name,closest.marker,Beta,LOD,FDR,has_interaction) 
C = C %>% dplyr::select(transcript,gene_name,closest.marker,Beta,LOD,FDR, has_interaction) 
colnames(B) = col_names_out
colnames(C) = col_names_out
out_list = list(A=one_pot_bulk_out, B=B,C=C)

openxlsx::write.xlsx(out_list,"tables/s7.xlsx")


