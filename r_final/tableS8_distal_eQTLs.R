A = cross_data$A$trans$hotspot_peaks_new
B = cross_data$B$trans$hotspot_peaks_new #$cis_test_with_disp_expr
C = cross_data$`3004`$trans$hotspot_peaks_new
A = A %>% filter(FDR < 0.05) %>% filter(chrom != tchr)
B = B %>% filter(FDR < 0.05) %>% filter(chrom != tchr)
C = C %>%  filter(FDR < 0.05) %>% filter(chrom != tchr)
A = A %>% inner_join(genes_name_trans,by=c("transcript"="gene_id"))
B = B %>% inner_join(genes_name_trans,by=c("transcript"="gene_id"))
C = C %>% inner_join(genes_name_trans,by=c("transcript"="gene_id"))


A_out = A %>% dplyr::select(transcript,gene_name, peak.marker, tchr,tpos, Beta.y,LOD,FDR,in.hotspot,bin,has_interaction_trans)
#A_out$LOD = 10^(-A_out$LOD)
B_out = B %>% dplyr::select(transcript,gene_name, peak.marker, tchr,tpos, Beta.y,LOD,FDR,in.hotspot,bin,has_interaction_trans)
#B_out$LOD = 10^(-B_out$LOD)
C_out = C %>% dplyr::select(transcript,gene_name, peak.marker, tchr,tpos, Beta.y,LOD,FDR,in.hotspot,bin,has_interaction_trans)
#C_out$LOD = 10^(-C_out$LOD)
#col_names_out = c("transcript","gene name","closest marker","estimate (one-pot)","p-value (one-pot)","adjusted p-value (one-pot)","has cell-cycle interaction")


col_names_table = c("transcript","gene name","peak marker","transcript chromosome","transcript position","estimate (one-pot)","LOD (one-pot)","FDR (one-pot)","eQTL in hotspot","hotspot bin","has cell-cycle interaction")
#hist(A_out$LOD)
colnames(A_out) = col_names_table
colnames(B_out) = col_names_table
colnames(C_out) = col_names_table
out_l = list(A=A_out,B=B_out,C=C_out)

openxlsx::write.xlsx(out_l,"tables/s8.xlsx")
