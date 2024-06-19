norm_diploid = function(cc_diploid){
  cc_diploid = cc_diploid %>% mutate(contrast2=str_replace_all(contrast," ","")) %>% mutate(contrast2=str_replace_all(contrast2,"\\(","")) %>%
    mutate(contrast2=str_replace_all(contrast2,"\\)","")) %>% mutate(contrast2=str_replace_all(contrast2,"/",":")) %>% 
    mutate(contrast2=str_replace_all(contrast2,"-"," ")) 
  
  cc_diploid = cc_diploid %>% mutate(contrast3=case_when(
    contrast2 == "G2:M S" ~ "S G2:M",
    contrast2 == "M:G1 S" ~ "S M:G1",
    TRUE ~ contrast2
  )) %>% mutate(estimate2=case_when(
    contrast2 == "G2:M S" ~  -estimate,
    contrast2 == "M:G1 S" ~ -estimate,
    TRUE ~ estimate
  ))
  
  cc_diploid = cc_diploid %>% mutate(p_adj = p.adjust(p.value,method="fdr"))
  return(cc_diploid)
}

load_cis_one_pot= function(cross){
  #sc_eqtl_row = summary_table %>% filter(cross == !!cross) %>% filter(type== "sceQTL")
  #print(sc_eqtl_row)
  cross_data_l = list()
  # Seurat objects
  cross_data_l[["one_pot_seurat"]] = list()
  for(ee in sets[[cross]]){
    experiment = experiments[ee]
    cc_in_f = paste(cc_dir,"/",experiment,"/cell_cycle_final.RDS",sep="")
    all_in_f = paste(cc_dir,"/",experiment,"/all_seurat.RDS",sep="")
    cc_in = readRDS(cc_in_f)
    all_in = readRDS(all_in_f)
    cross_data_l[["one_pot_seurat"]][[experiment]] = list(cc_in=cc_in,all_in=all_in)
  }
  # Get mapping results for cis # 
  in_combined_dir = paste0(combined_dir,"/",cross,"/")
  segdata = readRDS(paste0(in_combined_dir,"segData.RDS")) 
  chrom_pos = names(colMeans(segdata$Gsub))
  df_af = data.frame(chrom_pos=chrom_pos,means=colMeans(segdata$Gsub))
  #Chromosome III 198671..201177
 # print("GOT HER")
  #colM(segdata$Gsub)
  cross_data_l[["segData"]] = segdata 
  df_af = df_af %>% mutate(ss=str_split(chrom_pos,"_"))
  df_af = df_af %>% mutate(chrom=unlist(lapply(ss, function(x){x[1]})))
  df_af = df_af %>% mutate(pos=unlist(lapply(ss, function(x){x[2]})))
  df_af$pos = as.numeric(df_af$pos)
  cis_test_only_combined = readRDS(glue("{in_combined_dir}/cis_only_test_CombinedResults.RDS"))
  cis_test_interaction_betas = bind_rows(cis_test_only_combined,.id="cell_cycle")
  cis_test_interaction_betas$cell_cycle2=(str_replace_all(cis_test_interaction_betas$cell_cycle,":","/"))
  cis_test_only_combined_df = bind_rows(cis_test_only_combined,.id="cc")
  cis_test_only_combined_df = cis_test_only_combined_df %>% filter(cc != "combined") %>% group_by(transcript) %>% dplyr::summarise(Beta=get_theta(Beta,SE),SE=get_theta(SE,SE)) %>%
    inner_join(cis_test_only_combined$combined,by="transcript")
  sc_expr_df = data.frame(gene=names(colSums(segdata$Yr)),expr=colSums(segdata$Yr))
  sc_expr_df$cpm = sc_expr_df$expr/sum(sc_expr_df$expr)*1e6
  dispersion_seg = readRDS(glue("{in_combined_dir}/dispersions.RDS"))
  cis_test_with_disp_expr = cis_test_only_combined_df %>% inner_join(sc_expr_df,by=c("transcript"="gene")) %>% inner_join(dispersion_seg,by=c("transcript"="gene"))
  cut_cis = cut(log(abs(cis_test_with_disp_expr$Beta)), breaks=5)
  cis_test_with_disp_expr$cut_cis = cut_cis
  

  cis_test_only_combined_contrasts = readRDS(glue("{in_combined_dir}/cis_only_test_CombinedResultsContrasts.RDS"))
  contrasts_df = bind_rows(cis_test_only_combined_contrasts,.id="cc") %>% mutate(p_adj=p.adjust(pval,method="fdr"))
  cross_data_l[["cis_test_combined_contrasts"]] = cis_test_interaction_betas
  cross_data_l[["cis_test_interaction_betas"]] = contrasts_df
 
  sig_cis= contrasts_df %>% group_by(transcript)  %>% summarise(sig_count=sum(p_adj < 0.05))
  #p.value.cond
  
  cis_test_with_disp_expr$has_interaction = cis_test_with_disp_expr$transcript %in% sig_cis$transcript[sig_cis$sig_count > 0 ] 
  
  cross_data_l[["cis_test_with_disp_expr"]] = cis_test_with_disp_expr
  #cross_data_l[["hotspot_peaks"]] = readRDS(glue("{in_combined_dir}/hotspot_peaks.RDS"))
  cross_data_l[["fdr_cis_fx"]] = readRDS(glue("{in_combined_dir}/fdrfx_NB_combined.RDS"))
  
  sc_eqtl_row = summary_table %>% filter(cross == !!cross)  %>% filter(type!="eQTL")
  j = 1
  ase=list()
  for(in_folder_ase in sc_eqtl_row$folder){
    #print("HERE")
    #print(combined_dir)
    new_ase_folder = glue("{root_dir}/diploid_flip/")
    
    #new_ase_folder = "/media/theboocock/Data/Dropbox/PHDTHESIS/projects/single_cell_2021/rproj/out/diploid_flip/"
    old_ase_folder = glue("{root_dir}/") 
    #old_ase_folder = "/media/theboocock/Data/Dropbox/PHDTHESIS/projects/single_cell_2021/rproj/out/"
    #print(in_folder_ase)
    
    b_f = basename(in_folder_ase)
    
    
    
    in_folder_ase = glue("{new_ase_folder}{b_f}")
    #in_folder_cc = readRDS(glue("../rproj/out/cell_cycle/{b_f}"))
    cc_in = readRDS(glue("data/out/cell_cycle/{b_f}/{cross}/cell_cycle_final.RDS"))
    seurat_in = readRDS(glue("data//out/cell_cycle/{b_f}/{cross}/all.RDS"))
    type=sc_eqtl_row$type[j]  
    in_rds = readRDS(glue("{in_folder_ase}/bbin_{cross}.RDS"))
    #bbin_A_CCmanual.RDS
    in_cc_int = readRDS(glue("{in_folder_ase}/bbin_{cross}_CCmanual.RDS"))
    #readRDS(glue("{nbin_A.RDS")
    pair_df = attr(in_cc_int,"pairs") %>% bind_rows(.id="gene")
    pair_df = norm_diploid(pair_df)
    contrast_df=attr(in_cc_int,"contrasts") %>% bind_rows(.id="gene")
    in_rds$p_adj = p.adjust(in_rds$p.value,method = "fdr")
    #in_rds
    p2 = pair_df %>% filter(p_adj < 0.05) 
    with_int = unique(p2$gene)
    in_rds$has_interaction = in_rds$gene %in% with_int
    #ase_noise
    #in_cc_int %>% mutate()
    old_ase_folder= glue("{old_ase_folder}{b_f}") 
    
    #list.files(old_ase_folder)
    ase_noise_nbin1 = readRDS(glue("{old_ase_folder}/nbin1_{cross}.RDS"))
    
    ase_noise= readRDS(glue("{in_folder_ase}/{cross}-all_models.RDS"))
    ase_noise$has_interaction = ase_noise$gene %in% with_int
    nbin= readRDS(glue("{in_folder_ase}/nbin_{cross}.RDS"))
    gm = attr(nbin,'gm')
    gd = attr(nbin,"gd")
    gm_d = gm %>% inner_join(gd,by=c("gene","geno"),suffix=c(".mean",".disp")) 
    ase_noise = ase_noise %>% mutate(p_adj_disp = p.adjust(p.value.disp,"fdr"))
    ase_noise = ase_noise %>% mutate(p_adj_ase = p.adjust(p.value.cond,"fdr"))
    ase_data = readRDS(glue("{in_folder_ase}/aseData.RDS"))
    diploid_assignments = readRDS(glue("{in_folder_ase}/diploid_assignments.RDS"))
    phased_counts = readRDS(glue("{in_folder_ase}/phasedCounts_{cross}.RDS"))
    
    
    
    #if(cross %in% c("3004","B")){
    #  in_rds = in_rds %>% mutate(estimate=-estimate)
    #  in_rds = in_rds %>% mutate(statistic=-statistic)
    #  pair_df = pair_df   %>% mutate(estimate=-estimate, z.ratio=-z.ratio,estimate2=-estimate2)
    #  tmp = contrast_df$asymp.LCL
    #  contrast_df = contrast_df %>% mutate(emmean=-emmean,asymp.LCL=-asymp.UCL,asymp.UCL=-tmp)
    #  ase_noise %>% mutate(estimate=-estimate,statistic=-statistic,estimate.cond=-estimate.cond,
    #                       statistic.cond=-statistic.cond,estimate.disp=-estimate.disp,statistic.disp=-statistic.disp,
    #                       estimate.red.cond)
    #}
    ase_list = list(ase_cis=in_rds, ase_int=pair_df,ase_int_contrast=contrast_df, ase_noise=ase_noise,
                    ase_data=ase_data,diploid_assignments=diploid_assignments,
                    phased_counts=phased_counts,cc_in=cc_in,all_rds=seurat_in,geno_mean_disp=gm_d,ase_noise_nbin1=ase_noise_nbin1)
    
    if (type == "ASE"){
      cross_data_l[["ASE"]] = ase_list  
      # cross_data[[cross]][[""]]
    }else{
      cross_data_l[["ASE_REP"]] = ase_list  
    }
    j = j + 1
  }
  return(cross_data_l)
}
