

ap_icc  = readRDS("../rproj/out/combined/Ap/scICC.RDS")
genes = rownames(ap_icc)
ap_icc = ap_icc %>% as_data_frame()
ap_icc$gene_id = genes
ap_icc2 = ap_icc %>% inner_join(genes_name_trans,by=c("gene_id"))
colnames(ap_icc2) = c("ICC (cell-cycle)","ICC (broad-sense heritability H^2)","Raw cell-cycle variance","Raw genetic variance","Cell-cycle q-value","Genetic q-value","total umis","gene id","gene name")
ap_icc2 %>% write_tsv(file="tables/s5.tsv")
