#ASE_d = combined_objects$noise$ASE
#ASE_d %>% filter(gene, estimate, p.value, )
combined_objects$cis_eqtl_ase_only_with_onepot$sig = (combined_objects$cis_eqtl_ase_only_with_onepot$p_adj < 0.05)
#combined_objects$cis_eqtl_ase_only_with_onepot %>% group_by(cross) %>% summarise(sig=sum(sig), n=n())


out_ase = combined_objects$cis_eqtl_ase_only_with_onepot %>% inner_join(genes_name_trans,by=c("gene"="gene_id"))

out_ase = out_ase%>% dplyr::select(gene, gene_name,estimate,p.value,p_adj, has_interaction.x,Beta,LOD,FDR)#,cross)

sum(out_ase$has_interaction.x)
colnames(out_ase)
out_ase_col = c("transcript","gene name","estimate (ASE)","p-value (ASE)","adjust p-value (ASE)","has cell-cycle interaction (ASE)","estimate (one-pot)","LOD (one-pot)","FDR (one-pot)")
colnames(out_ase) = out_ase_col

l_ase = split(out_ase,f=combined_objects$cis_eqtl_ase_only_with_onepot$cross)
names(l_ase) = c("C","A","B")
l_ase = l_ase[c("A","B","C")]
openxlsx::write.xlsx(l_ase,file="tables/s9.xlsx")

#nn#ASE_d %>% filter(cross ==)
