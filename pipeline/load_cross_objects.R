source("pipeline/packages.R")
source("pipeline/qtl_utils.R")
source("pipeline/enrichment_utils.R")
source("pipeline/annotation_utils.R")
source("pipeline/plot_fx.R")
cc_dir = normalizePath("data//out/cell_cycle/")
combined_dir = normalizePath("data//out/combined/")


summary_table = read_tsv("tables/s4.csv")
dir.create("figure_drafts")
fig.dir = "figure_drafts"
cross_data= list()
#i = 1
plot=F

source("pipeline/cis_objects.R")
source("pipeline/trans_objects.R")
source("pipeline/vars.R")

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

#A (BY-, RM+) # TODO: 
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
source("load_single_cell_gpa1.R")
# TODO don't load things that take ages to compute just load the objects and break them out.
source("load_gpa1_phenotype_and_mating_data.R")
source("load_gpa1_phylogenetic_data.R")

combined_objects = list()
combined_objects[["noise"]] = list()

#hotspot_peaksA = cross_data$A$trans$combined_peaks %>% mutate(cross="A")
#hotspot_peaks3004 = cross_data$`3004`$trans$combined_peaks %>% mutate(cross="3004")
#hotspot_peaksB = cross_data$B$trans$combined_peaks %>% mutate(cross ="B")

for(cross in c("A","B","3004")){
  
  ## ZIP up the region_df file ## can do it here for now.
  print(cross)
  cis = cross_data[[cross]]$cis$cis_test_with_disp_expr %>% mutate(cross = !!cross)
  trans_peaks = cross_data[[cross]]$trans$hotspot_peaks_new %>% mutate(cross = !!cross)
  
  cis_ase_only = cross_data[[cross]]$cis$ASE$ase_cis %>% mutate(cross = !!cross)
  
  cis_ase = cross_data[[cross]]$cis$ASE$ase_noise %>% mutate(cross = !!cross)
  #cis_ase_only_m
  #cor(cis_ase$estimate,cis_ase_only$estimate)
  cis_ase_only_m = left_join(cis_ase_only,cis,by=c("gene"="transcript","cross"))
  cis_noise_m  = inner_join(cis,cis_ase,by=c("transcript"="gene","cross")) 
  print(nrow(cis_ase_only_m))
  ### Probably flip these ###
 # if (cross %in% c("3004","B")){
  #  cis_noise_m$Beta = -cis_noise_m$Beta
    #cis_ase_only$estimate = -cis_ase_only$estimate
  #}
  
  
  ccl = cross_data[[cross]]$trans$cell_cycle_lods %>% as_data_frame()
  combined_objects[["cell_cycle_lods"]] = rbind(combined_objects[["cell_cycle_lods"]],ccl)
  #combined_objects[["cis_ase"]]
  combined_objects[["cis_eqtl_ase_only_with_onepot"]] = rbind(combined_objects[["cis_eqtl_ase_only_with_onepot"]],cis_ase_only_m)
  combined_objects[["cis_eqtl_ase_noise_with_onepot"]] = rbind(combined_objects[["cis_eqtl_ase_noise_with_onepot"]],cis_noise_m)
  
  combined_objects[["cis_table"]] = rbind(combined_objects[["cis_table"]],cis)
  combined_objects[["hotspot_peaks"]] = rbind(combined_objects[["hotspot_peaks"]],trans_peaks)
  combined_objects[["cis_eqtl_ase"]] = rbind(combined_objects[["cis_eqtl_ase"]],cis_ase_only)
  
  
  n1 = cross_data[[cross]]$cis$ASE$ase_noise %>% mutate(cross = !!cross)
  #n2 = cross_data[[cross]]$cis$ASE_REP$ase_noise %>% mutate(cross = !!cross)
  
  #n1_rep_m = (n1 %>% left_join(n2,by=c("gene","cross"),suffix=c(".ASE",".REP")))
  merged_with_cis_noise = cis %>% full_join(n1,by=c("transcript"="gene"),suffix=c(".CIS",".ASE"))
  
  #print(table(merged_with_cis_noise$cross.ASE))
  combined_objects[["noise"]][["NOISE_CIS_M"]]=rbind(combined_objects[["noise"]][["NOISE_CIS_M"]],merged_with_cis_noise)
  combined_objects[["noise"]][["ASE"]] = rbind(combined_objects[["noise"]][["ASE"]],n1)
  #combined_objects[["noise"]][["ASE_REP"]] = rbind(combined_objects[["noise"]][["ASE_REP"]],n2)
  #combined_objects[["noise"]][["ASE_M"]]= rbind(combined_objects[["noise"]][["ASE_M"]],n1_rep_m)
  
  noise_gm = cross_data[[cross]]$cis$ASE$geno_mean_disp %>% mutate(cross = !!cross)
  
  combined_objects[["noise"]][["ASE_EMMEANS"]] = rbind(combined_objects[["noise"]][["ASE_EMMEANS"]],noise_gm)
  
  #ccl = cross_data[[cross]]$trans$cell_cycle_lods
  #ccl = ccl %>% mutate(cross = !!cross) %>% as_data_frame()
  #combined_objects[["cell_cycle_lods"]] = rbind(combined_objects[["cell_cycle_lods"]],ccl)
  #cross_data$A$cis$ASE$ase_int
  
  #cross_data$A$cis$ASE$geno_mean_disp
  
}

combined_objects$noise$ASE  = combined_objects$noise$ASE %>% mutate(sig_new_filt = ( combined_objects$noise$ASE$estimate.cond*combined_objects$noise$ASE$estimate.disp < 0) & p_adj_disp < .05)# %>%
combined_objects$noise$ASE_EMMEANS = combined_objects$noise$ASE_EMMEANS  %>% mutate(cross2  = case_when(
  cross == "3004" ~ "CBS2888 (red) x YJM981 (purple)",
  cross == "B" ~ "YJM145 (red) x YPS163 (purple)",
  cross == "A" ~ "BY (red) x RM (purple)"
))

combined_objects$cis_eqtl_ase_noise_with_onepot = combined_objects$cis_eqtl_ase_noise_with_onepot %>% mutate(cross2  = case_when(
  cross == "3004" ~ "CBS2888 x YJM981 (C)",
  cross == "B" ~ "YJM145 x YPS163 (B)",
  cross == "A" ~ "BY x RM (A)"
))

combined_objects$cis_eqtl_ase_only_with_onepot = combined_objects$cis_eqtl_ase_only_with_onepot %>% mutate(cross2  = case_when(
  cross == "3004" ~ "CBS2888 x YJM981 (C)",
  cross == "B" ~ "YJM145 x YPS163 (B)",
  cross == "A" ~ "BY x RM (A)"
))

#combined_objects[["trans_maps"]] = rbind(hotspot_peaksA,hotspot_peaksB,hotspot_peaks3004)

#cis_A = cross_data$A$cis$cis_test_with_disp_expr %>% mutate(cross = "A")

cross_data$A$trans$hotspot_table$cross = "A"
cross_data$B$trans$hotspot_table$cross = "B"
cross_data$`3004`$trans$hotspot_table$cross = "3004"

A_hot = cross_data$A$trans$hotspot_bins_final #%>% filter(count > cross_data$A$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()
B_hot = cross_data$B$trans$hotspot_bins_final# %>% filter(count > cross_data$B$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()
C_hot = cross_data$`3004`$trans$hotspot_bins_final#
A_old_hot = cross_data$A_bulk$trans$hotspot_bins_final# %>% filter(count > cross_data$A_bulk$trans$hotspot_threshold) %>% mutate(seqnames=chr,start=pos,width=1) %>% as_granges()



