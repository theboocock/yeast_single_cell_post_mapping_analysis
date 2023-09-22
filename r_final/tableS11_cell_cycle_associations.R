cc_lods2 = combined_objects$cell_cycle_lods %>% dplyr::select(seqnames,marker, ci_left, ci_right, beta,lod)
cc_lod_str = c("chromosome","C.I. left","C.I. right","estimate (cell-cycle effect)","LOD (cell-cycle effect)")
c_l = split(cc_lods2,combined_objects$cell_cycle_lods$experiment)
names(c_l) = c("C","A","B")
c_l = c_l[c("A","B","C")]
openxlsx::write.xlsx(c_l, file="tables/s11.xlsx")
