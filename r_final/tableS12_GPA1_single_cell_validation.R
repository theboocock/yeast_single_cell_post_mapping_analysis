#bl
# A  =420/421 82R
# B = 416/417 WT
pm = cross_data$B$trans$hotspot_enrichments_and_overlaps$`chrVIII:46887-140660_51`$directional_hotspot_distal
#pm$gene_name[!]
blm=mk%>% rownames_to_column(var="gene")
bl = pm %>% left_join(blm,by=c("gene_name"="gene"))

bl2 = bl %>% dplyr::select(transcript, gene_name,peak.marker,Beta.y, LOD,FDR,avg_log2FC,p_val,p_val_adj) %>% mutate(avg_log2FC=-avg_log2FC)
sum(bl2$Beta.y * bl2$avg_log2FC > 0)

colnames(bl2) = c("transcript","gene name","peak marker","estimate (one-pot)", "LOD (one-pot)","FDR (one-pot)","estimate (single-cell validation)","p-value (single-cell validation)","adjusted p-value (single-cell validation)")
openxlsx::write.xlsx(bl2, file="tables/s12.xlsx")
