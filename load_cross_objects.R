source("packages.R")
source("qtl_utils.R")
source("enrichment_utils.R")
source("annotation_utils.R")

cc_dir = normalizePath("../rproj/out/cell_cycle/")
combined_dir = normalizePath("../rproj/out/combined/")


summary_table = read_csv("../figure_and_tables_paper/Tables/s4.csv")
dir.create("figure_drafts")
fig.dir = "figure_drafts"
cross_data= list()
#i = 1
plot=F

source("cis_objects.R")
source("trans_objects.R")
source("vars.R")

#library(foreach)
#cl <- parallel::makeForkCluster(4)
#doParallel::registerDoParallel(cl)

#ret_cross =foreach(cross = unique(summary_table$cross),.combine="c") %dopar% {
  #print(cross)
  #cross_tmp = list()
  #if(cross == "Ap"){
  #  return(cross_tmp)
  #}
#  if(cross != "Ap"){
#    cross_tmp[["cis"]] = load_cis_one_pot(cross)
#    cross_tmp[["trans"]]= load_trans_one_pot(cross)
    
    # Get HAPLOIDS
    ###
#    if(cross == "A"){
#      ## Do the comparsions ##
#    }
    
    
    # GET DIPLOIDS
    
    
#  }
#  return(cross_tmp)
#}

#A (BY-, RM+)
#3004 (CBS-, YJM981+)
#B (YJM145-, YPS163+

for(cross in unique(summary_table$cross)){
  print(cross)
  if(cross == "Ap"){
    next
  }
  if(cross != "Ap"){
    cross_data[[cross]][["cis"]] = load_cis_one_pot(cross)
    cross_data[[cross]][["trans"]]= load_trans_one_pot(cross)
    
    
    # Get HAPLOIDS
    ###
    if(cross == "A"){
      ## Do the comparsions ##
    }
    
    
    # GET DIPLOIDS
    
    
  }
  #break
  #i = 1
}
##### GET THE BULK STUFF ######
source("load_by_rm_bulk.R")

combined_objects = list()
combined_objects[["noise"]] = list()

#hotspot_peaksA = cross_data$A$trans$combined_peaks %>% mutate(cross="A")
#hotspot_peaks3004 = cross_data$`3004`$trans$combined_peaks %>% mutate(cross="3004")
#hotspot_peaksB = cross_data$B$trans$combined_peaks %>% mutate(cross ="B")

for(cross in c("A","B","3004")){
  print(cross)
  cis = cross_data[[cross]]$cis$cis_test_with_disp_expr %>% mutate(cross = !!cross)
  trans_peaks = cross_data[[cross]]$trans$hotspot_peaks_new %>% mutate(cross = !!cross)
  
  cis_ase_only = cross_data[[cross]]$cis$ASE$ase_cis %>% mutate(cross = !!cross)
  
  cis_ase = cross_data[[cross]]$cis$ASE$ase_noise %>% mutate(cross = !!cross)
  
  
  cis_noise_m  = inner_join(cis,cis_ase,by=c("transcript"="gene","cross")) 
  if (cross %in% c("3004","B")){
    cis_noise_m$Beta = -cis_noise_m$Beta
    #cis_ase_only$estimate = -cis_ase_only$estimate
  }
  #combined_objects[["cis_ase"]]
  combined_objects[["cis_eqtl_ase_noise"]] = rbind(combined_objects[["cis_eqtl_ase_noise"]],cis_noise_m)
  
  combined_objects[["cis_table"]] = rbind(combined_objects[["cis_table"]],cis)
  combined_objects[["hotspot_peaks"]] = rbind(combined_objects[["hotspot_peaks"]],trans_peaks)
  combined_objects[["cis_eqtl_ase"]] = rbind(combined_objects[["cis_eqtl_ase"]],cis_ase_only)
  
  
  n1 = cross_data[[cross]]$cis$ASE$ase_noise %>% mutate(cross = !!cross)
  n2 = cross_data[[cross]]$cis$ASE_REP$ase_noise %>% mutate(cross = !!cross)
  
  n1_rep_m = (n1 %>% left_join(n2,by=c("gene","cross"),suffix=c(".ASE",".REP")))
  merged_with_cis_noise = cis %>% full_join(n1,by=c("transcript"="gene"),suffix=c(".CIS",".ASE"))
  
  #print(table(merged_with_cis_noise$cross.ASE))
  combined_objects[["noise"]][["NOISE_CIS_M"]]=rbind(combined_objects[["noise"]][["NOISE_CIS_M"]],merged_with_cis_noise)
  
  combined_objects[["noise"]][["ASE"]] = rbind(combined_objects[["noise"]][["ASE"]],n1)
  combined_objects[["noise"]][["ASE_REP"]] = rbind(combined_objects[["noise"]][["ASE_REP"]],n2)
  combined_objects[["noise"]][["ASE_M"]]= rbind(combined_objects[["noise"]][["ASE_M"]],n1_rep_m)
  
  noise_gm = cross_data[[cross]]$cis$ASE$geno_mean_disp %>% mutate(cross = !!cross)
  
  combined_objects[["noise"]][["ASE_EMMEANS"]] = rbind(combined_objects[["noise"]][["ASE_EMMEANS"]],noise_gm)
  
  ccl = cross_data[[cross]]$trans$cell_cycle_lods
  ccl = ccl %>% mutate(cross = !!cross) %>% as_data_frame()
  combined_objects[["cell_cycle_lods"]] = rbind(combined_objects[["cell_cycle_lods"]],ccl)
  #cross_data$A$cis$ASE$ase_int
  
  #cross_data$A$cis$ASE$geno_mean_disp
  
}

combined_objects$noise$ASE  = combined_objects$noise$ASE %>% mutate(sig_new_filt = (p_adj_ase > .05 | (combined_objects$noise$ASE$estimate.cond*combined_objects$noise$ASE$estimate.disp < 0)) & p_adj_disp < 1e-4)# %>%
combined_objects$noise$ASE_EMMEANS = combined_objects$noise$ASE_EMMEANS  %>% mutate(cross2  = case_when(
  cross == "3004" ~ "YJM981 (red) x CBS2888 (purple)",
  cross == "B" ~ "YPS163 (red) x YJM145 (purple)",
  cross == "A" ~ "BY (red) x RM (purple)"
))

combined_objects$cis_eqtl_ase_noise = combined_objects$cis_eqtl_ase_noise %>% mutate(cross2  = case_when(
  cross == "3004" ~ "YJM981 x CBS2888",
  cross == "B" ~ "YPS163 x YJM145",
  cross == "A" ~ "BY x RM"
))

#combined_objects[["trans_maps"]] = rbind(hotspot_peaksA,hotspot_peaksB,hotspot_peaks3004)

#cis_A = cross_data$A$cis$cis_test_with_disp_expr %>% mutate(cross = "A")





