process_hotspot_and_combined = function(hotspot_peaks, combined_peaks){
  hotspot_peaks= bind_rows(hotspot_peaks,.id = "cell_cycle_stage")
  hotspot_peaks= hotspot_peaks %>% dplyr::filter(cell_cycle_stage == "combined")
  str_split_bin =  str_split(hotspot_peaks$bin,":")
  pos_string = unlist(lapply(str_split_bin, function(x){x[[2]]}))
  pos_string_l = unlist(lapply(str_split(pos_string,"-"), function(x){x[1]}))
  hotspot_peaks$bin_pos_l = factor(as.numeric(pos_string_l))
  
  combined_peaks = combined_peaks %>% mutate(
    Z_G1 = Beta_G1/SE_G1,
    `Z_G2:M` = `Beta_G2:M`/`SE_G2:M`,
    `Z_M:G1` = `Beta_M:G1`/`SE_M:G1`,
    `Z_G1:S` = `Beta_G1:S`/`SE_G1:S`,
    `Z_S` = Beta_S /SE_S
  )
  # Add betas to hotspot peaks #
  x_beta= hotspot_peaks %>% pivot_longer(cols=starts_with("Beta_"),values_to = "beta")# %>% group_by(transcript,peak.marker) %>% summarise(Beta=get_theta(value))
  x_se = hotspot_peaks %>% pivot_longer(cols=starts_with("SE_", ), values_to = "se")
  x_beta$se = x_se$se
  x_b2 = x_beta %>% group_by(chrom,transcript,peak.marker) %>% summarise(Beta=get_theta(beta,se),SE=get_theta(se,se))
  #hotspot_peaks$Beta = x_b2$Beta
  #hotspot_peaks$SE = x_b2$SE 
  hotspot_peaks = hotspot_peaks %>% inner_join(x_b2,by=c("chrom","transcript","peak.marker"))
  hotspot_peaks$chrom_short = str_replace_all(hotspot_peaks$chr,"chr","")
  hotspot_peaks$tchrom_short = str_replace_all(hotspot_peaks$tchr,"chr","")
  combined_peaks$chrom_short = str_replace_all(combined_peaks$chr,"chr","")
  combined_peaks$tchrom_short = str_replace_all(combined_peaks$tchr,"chr","")
  
  order_levels = c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")
  rev_order_levels = rev(order_levels)
  hotspot_peaks$chrom_short_f = factor(hotspot_peaks$chrom_short,levels=order_levels)
  combined_peaks$chrom_short_f = factor(combined_peaks$chrom_short,levels=order_levels)
  hotspot_peaks$tchrom_short_f = factor(hotspot_peaks$tchrom_short,levels=rev_order_levels)
  combined_peaks$tchrom_short_f = factor(combined_peaks$tchrom_short,levels=rev_order_levels)
  
  return(list(combined_peaks=combined_peaks,hotspot_peaks=hotspot_peaks))
}

#"get_lod_from_r"


load_trans_one_pot = function(cross){
  sc_eqtl_row = summary_table %>% filter(cross == !!cross) %>% filter(type== "sceQTL")
  cross_data_l = list()
  in_combined_dir = paste0(combined_dir,"/",cross,"/")
  hotspot_peaks = readRDS(glue("{in_combined_dir}/hotspot_peaks.RDS"))
  combined_peaks = readRDS(glue("{in_combined_dir}/LOD_NB_combined_peaks.RDS"))
  segdata = readRDS(paste0(in_combined_dir,"segData.RDS")) 
  cross_data_l[["segdata"]] = segdata
  hotspot_threshold = attr(hotspot_peaks$combined,"threshold")
  ret = process_hotspot_and_combined(hotspot_peaks,combined_peaks)
  cross_data_l[["combined_peaks"]] =ret$combined_peaks
 # cross_data_l[["hotspot_peaks"]] = ret$hotspot_peaks
  cross_data_l[["hotspot_threshold"]] = hotspot_threshold
  cc_in = glue("{in_combined_dir}/cell_cycle_assignment_LOD.RDS")
  cclod = readRDS(cc_in)
  fwer_threshcc = attr(cclod,"FWER_thresh_5%")
  cc_b = glue("{in_combined_dir}/cell_cycle_assignment_Betas.RDS")
  betas = readRDS(cc_b)
  ses = readRDS(glue("{in_combined_dir}/cell_cycle_assignment_SEs.RDS"))
  cc_row_names = rownames(betas)
  cell_cycle_rds_df = make_cc_lod_df(cclod,experiment = cross)
  lod_df = lod_drop_cc(cell_cycle_rds_df,lod_threshold = fwer_threshcc,lod_drop = 1.5)
  cell_cycle_beta_df = make_cc_lod_df(betas,experiment = cross) %>% dplyr::rename(beta=lod)
  cell_cycle_se_df = make_cc_lod_df(ses,experiment = cross) %>% dplyr::rename(se=lod)
  cell_cycle_beta_df2 = cell_cycle_beta_df %>% mutate(se=cell_cycle_se_df$se) # %>% select(cell_cycle,marker,beta,se
  cell_cycle_beta_df3 = cell_cycle_beta_df2 %>% dplyr::select(cell_cycle,marker,beta,se)
  lod_df_beta = lod_df %>% inner_join(cell_cycle_beta_df3,by=c("cell_cycle","marker")) 
  lod_df_beta = lod_df_beta %>% mutate(ci_left=ifelse(is.na(ci_left ),1,ci_left))
  #is.na()
  chrom_tmp = lod_df_beta$chrom[which(is.na(lod_df_beta$ci_right))]
  lengths = chrom_lengths$lengths[match(chrom_tmp,chrom_lengths$chrom)]
  lod_df_beta$ci_right[is.na(lod_df_beta$ci_right)] = lengths
  lod_df_beta = lod_df_beta %>% filter(chrom != "chrIII")
  #lod_df_beta = lod_df_beta %>% mutate(ci_left=ifelse(is.na(ci_left ),1,ci_left))
  lod_df_beta = lod_df_beta %>% mutate(seqnames=chrom, start=ci_left, end=ci_right) %>% as_granges()
  lod_df2 = lod_df_beta %>% reduce_ranges() #%>% 
  ab = findOverlaps(lod_df_beta,lod_df2)  
  bin_df = lod_df2[subjectHits(ab)] %>% as_data_frame() %>% mutate(bin=paste(seqnames,":",start,"-",end,sep=""))
  lod_df_beta$bin = bin_df$bin
  
  cell_cycle_beta_df2 = cell_cycle_beta_df2 %>% inner_join(cell_cycle_rds_df %>% dplyr::select(cell_cycle,marker,lod),by=c("cell_cycle","marker"))
  cell_cycle_beta_df2$chrom_short_f = convert_chrom_to_simple_factor(cell_cycle_beta_df2$chrom)
  
  cross_data_l[["cell_cycle_beta_all"]] = cell_cycle_beta_df2
  
  cross_data_l[["cell_cycle_lods"]] = lod_df_beta
  cross_data_l[["cell_cycle_threshold"]] = fwer_threshcc
  ### Get bins etc for each of the crosses to match what josh generated for plotting ###
  
  cbin50k = getcbin50k_genome()
  
  #cbin50k = getcbin50k(segdata$Gsub)
  hotspot_table = makeBinTable_hotspot(ret$hotspot_peaks,cbin50k)
  sig.hp=qpois(1-(.05/length(hotspot_table$pos)),ceiling(mean(hotspot_table$count)))+1
  markerGR = getMarkerGRanges(list(t(segdata$Gsub)))
  hotspot_peaks_n=getBinAssignment(ret$hotspot_peaks, cbin50k,markerGR )
  hotspot_peaks_n$in.hotspot = F
  sig.hp.names = table(hotspot_peaks_n$bin)  > sig.hp
  hotspot_peaks_n$in.hotspot = sig.hp.names[hotspot_peaks_n$bin]
  cross_data_l[["hotspot_threshold"]] = sig.hp
  
  ### TODO add cell-cycle interactions ####
  cross_data_l[["hotspot_table"]] = hotspot_table
  cross_data_l[["hotspot_table"]]$chrom_f = convert_chrom_to_simple_factor(cross_data_l[["hotspot_table"]]$chr)
  cell_cycle_interactions = readRDS(glue("{in_combined_dir}/trans_CC_test_CombinedResultsContrasts.RDS"))
  trans_cc = bind_rows(cell_cycle_interactions,.id="cc")
  cross_data_l[["fdr_cis_fx"]] = readRDS(glue("{in_combined_dir}/fdrfx_NB_combined.RDS"))
  
  trans_cc_df = trans_cc %>% mutate(p_adj =p.adjust(pval,method="fdr"))
  #trans_cc_df %>% filter(p_adj < 0.05)
  trans_cc_df_wider = trans_cc_df %>% pivot_wider(id_cols = c("chrom","peak.marker","transcript"),names_from = c("cc"),values_from = c("Zstat","p_adj"))
  cross_data_l[["trans_cc_df"]] = trans_cc_df
  cross_data_l[["trans_cc_df_wider"]] = trans_cc_df_wider
  #ret$hotspot_peaks
  cross_data_l[["trans_cc_df"]]$chrom_f = convert_chrom_to_simple_factor(cross_data_l[["trans_cc_df"]]$chrom)
  cross_data_l[["trans_cc_df_wider"]]$chrom_f = convert_chrom_to_simple_factor(cross_data_l[["trans_cc_df_wider"]]$chrom)
  
  #aa$bin  == ret$hotspot_peaks$bin
  
  #### add chromosome simple ##
  
  
  with_interactions=trans_cc_df %>% group_by(transcript,peak.marker) %>% summarise(p_point_o_five=sum(p_adj < 0.05))
  with_interactions = with_interactions %>% mutate(has_interaction_trans=p_point_o_five > 0)
  #with_interactions
  #hotspot_pe
  hotspot_peaks_n = hotspot_peaks_n %>% inner_join(with_interactions,by=c("transcript","peak.marker"))
 # has_interaction_trans = hotspot_peaks_n$transcript %in% with_interactions$transcript[with_interactions$p_point_o_five >0]
  #  ret$hotspot_peaks$has_interaction = has_interaction_trans
  #hotspot_peaks_n$has_interaction = has_interaction_trans
  
  cross_data_l[["hotspot_peaks_old"]] = ret$hotspot_peaks
  cross_data_l[["hotspot_peaks_new"]] = hotspot_peaks_n 
  hotspot_table_cc = makeBinTable_hotspot_cc(hotspot_peaks_n,cbin50k)
  ##  TODO: Add the start and end of chromosomes .... for plottinng #
  ### Identify hotspots ###
  fdr_fx = readRDS(glue("{in_combined_dir}/fdrfx_NB_combined.RDS"))
  #hotspot_table_cc  %>% filter(int_count > fwer_threshcc) 
  background = colnames(segdata$Yr)
  
  #cross = cross
  #hotspot_peaks = hotspot_peaks_n
  #combined_peaks = ret$combined_peaks
  
    
  
  enrich_list = hotspot_enrichment_and_function(cross=cross,
                                                hotspot_peaks =  hotspot_peaks_n,
                                                combined_peaks = ret$combined_peaks,
                                                lod_df_beta = lod_df_beta,
                                                cell_cycle_beta_df2 = cell_cycle_beta_df2,
                                                segdata = segdata,
                                                fdr_fx = fdr_fx,
                                                background=background)
  cross_data_l[["hotspot_list"]] = enrich_list$hotspot_list
  
  
  
  hotspot_bins_final = enrich_list$hotspot_list %>% group_by(bin,hotspot_pos,chrom) %>% summarise(n=n())
  
  hotspot_bins_final$chrom = unlist(lapply(str_split(hotspot_bins_final$bin,":"), function(x){x[1]}))
  pos_str = unlist(lapply(str_split(hotspot_bins_final$bin,":"), function(x){x[2]}))
  start = as.numeric(unlist(lapply(str_split(pos_str,"-"), function(x){x[1]})))
  end =as.numeric(unlist(lapply(str_split(pos_str,"-"), function(x){x[2]})))
  hotspot_bins_final$start = start 
  hotspot_bins_final$end = end
  hotspot_bins_final_gr = makeGRangesFromDataFrame(hotspot_bins_final)
  hotspot_bins_final_gr$n = hotspot_bins_final$n
  hotspot_bins_final_gr$bin = hotspot_bins_final$bin
  hotspot_bins_final_gr$hotspot_pos = hotspot_bins_final$hotspot_pos
  cross_data_l[["hotspot_bins_final"]] = hotspot_bins_final_gr
  cross_data_l[["hotspot_list"]]$chrom_f = convert_chrom_to_simple_factor(cross_data_l[["hotspot_list"]]$chrom)
  
  cross_data_l[["hotspot_enrichments_and_overlaps"]] = enrich_list$annotation_list
  #ret$hotspot_peaks
  #hotspot_table_cc = makeBinTable_hotspot_cc(enrich_list$hotspot_list,cbin50k)
  
  in_base_dir="/media/theboocock/Data/Dropbox/PHDTHESIS/projects/single_cell_2021/rproj/out/cell_cycle/"
  
  cross_data_l[["seurat_objects_cc"]] = list()
  for(idxs in sets[[cross]]){
    print(idxs)
    in_cc_path = names(cList)[idxs]
    cell_cycle_rds = readRDS(paste(in_base_dir,names(cList)[idxs],"/cell_cycle_final.RDS",sep=""))
    #segdata$barcode.features
    bc = segdata$cell.cycle.df$cell_name[(segdata$cell.cycle.df$named_dataset == in_cc_path)]
    cell_cycle_rds$cc_seurat = cell_cycle_rds$cc_seurat[,bc]
    cross_data_l[["seurat_objects_cc"]][[in_cc_path]] = cell_cycle_rds
    #readRDS()
  }
  
  return(cross_data_l)
}

