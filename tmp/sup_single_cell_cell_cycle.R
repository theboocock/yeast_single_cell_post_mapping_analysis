sig_vc = readRDS("../rproj/out/combined/Ap/scICC_3VC.RDS")
#ap_vc = readRDS("../rproj/out/combined/Ap/scVCmodels.RDS")
scvc  =readRDS("../rproj/out/combined/Ap/scVCmodels.RDS")
ap = readRDS("../rproj/out/combined/Ap/scICC.RDS")
genes = rownames(ap)
ap_var = ap %>% as_data_frame()
ap_var$gene = genes 
fisher.test(table(ap_var$gene %in% cell_cycle_big_df$ORF,ap_var$CC.q < 0.05))
