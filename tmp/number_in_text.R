

### ###
summary_table
#

#  single-cell RNA-seq to obtain the transcriptomes of 7,124 single cells (Methods, Table S4). We captured a median of 1,514 unique RNA molecules (UMIs) and a median of 1,091 expressed SNPs per cell (Table S4).

summary_table[1,]
# We classified individual haploid yeast cells into five different stages of the 
#cell cycle (M/G1, G1, G1/S, S, G2/M) via unsupervised clustering of the expression of 787 cell-cycle-regulated genes (Spellman et al. 1998) 
#in combination with 22 cell-cycle-informative marker genes
#  Using this classification approach, we found that 2,139 genes displayed significant variation by cell cycle stage (likelihood ratio test, FDR < 0.05; Table S5
sum(ap_icc$CC.q < 0.05,na.rm=T)
# We observed a median of 17 cells per segregant, with 277 of the segregants sampled more than ten times (Figure S4)
median(seg_match_df$Freq)
sum(seg_match_df$Freq > 10)
#The genotypes measured from scRNA-seq data were in high agreement with those obtained from whole-genome sequencing of the same strains (median genotype agreement 92.5%). T
median(ap_df$barcode.features$accuracy)
#median(accuracy)

nrow(ap_combined %>% filter(p_adj.old < 0.05))
nrow(ap_combined %>% filter(p_adj.hmm < 0.05))
sum(ap_combined$sig_both)

#We sampled 250,000 MATa segregants and used scRNA-seq to measure the expression of 5,435 transcripts in 27,744 single cells
summary_table[2,]

## . We mapped 1,031 local eQTLs at a FDR of 5%. We compared the results in the BYxRM cross to those based on expression measurements from bulk RNA-seq (Albert et al. 2018) and found that 717 (69.5%) of the 1,031 local eQTLs
#were also detected as statistically significant in that study, and an additional 108 local eQTLs showed effects in the same direction (Figure 1C, 2A).

sum(one_pot_bulk$FDR < 0.05,na.rm=T)
sum(one_pot_bulk$FDR < 0.05 & one_pot_bulk$Beta * one_pot_bulk$coefficient > 0 &one_pot_bulk$p_adj < 0.05,na.rm=T )
717/1031
##
sum(one_pot_bulk$FDR < 0.05 & one_pot_bulk$Beta * one_pot_bulk$coefficient > 0,na.rm=T)
825/1031



###
### We mapped 1,562 distal eQTLs at a FDR of 5% (Figure 2A)(Table SX)
nrow(cross_data$A$trans$hotspot_peaks_new)
## We identified 12 trans eQTL hotspot loci in our one-pot eQTL experiment and bulk eQTL mapping 
##in segregants from the same cross identified 21 hotspot loci using the same criteria (Figure 2B). 
A_hot
A_old_hot

#### Five of the trans-eQTL hotspot loci identified by one-pot eQTL in segregants from 
##the BY and RM cross overlapped hotspots seen in the bulk eQTL mapping including the well described MKT1 
#(Zhu et al. 2008), GPA1 (Yvert et al. 2003), IRA2 (Smith and Kruglyak 2008), and HAP1 (Brem et al. 2002) hotspots.


m_hotspots= A_hot %>% join_overlap_left(A_old_hot)
m_hotspots_old = A_old_hot %>% join_overlap_left(A_hot)
m_hotspots_old[m_hotspots_old$hotspot_pos.x == 367656,]

m_hotspots[!is.na(m_hotspots$bin.y),]


m_hotspots[is.na(m_hotspots$bin.y),]
## Functional characterization of novel trans eQTL hotspots in the BY and RM cross ##
m_hotspots[is.na(m_hotspots$bin.y),]


## One-pot eQTL mapping in new crosses ###

#We took advantage of the convenience of one-pot eQTL mapping and applied it to two additional yeast crosses, 
#one between the strains YJM981 and CBS2888 (6,595 cells, 4,571 transcripts)
#, and another between the strains YPS163 and YJM145 (14,823 cells, 5,091 transcripts; Supplementary Note QC). 

summary_table[3,]
# TODO: fix with new table
summary_table[4,]

#cross_data$B$trans$combined_peaks %>% filter(FDR < 0.05)

#  We classified every cell from these crosses into the 5 different stages of the cell cycle as described above (Figure S15-S17). 
#We mapped 1,914 local eQTLs (721 local eQTLs for YJM981xCBS2888, and 1,193 local eQTLs for YPS163xYJM145)
nrow(cross_data$B$cis$cis_test_with_disp_expr %>% filter(FDR < 0.05) )
nrow(cross_data$`3004`$cis$cis_test_with_disp_expr %>% filter(FDR < 0.05))
# nd 1,193 local eQTLs for YPS163xYJM145) and 1,626 distal eQTLs (1,126 distal eQTLs for YJM981xCBS2888, 
#and 550 distal eQTLs for YPS163xYJM145) at a FDR of 5% (Figure 3A-B).
nrow(cross_data$B$trans$hotspot_peaks_new)
nrow(cross_data$`3004`$trans$hotspot_peaks_new)


## ##

B_hot
C_hot

## ##

## Noise eQTLs ###

sum(combined_objects$noise$ASE$p_adj_disp < 0.05,na.rm=T)
combined_objects$noise$ASE %>% group_by(cross) %>% summarise(sig=sum(p_adj_disp < 0.05,na.rm=T),n=n())

sum(combined_objects$noise$ASE$sig_new_filt,na.rm=T)

combined_objects$noise$ASE %>% group_by(cross) %>% summarise(sig=sum(sig_new_filt,na.rm=T),n=n())



#  We identified that 116 (3.9%) of our 2,945 local eQTLs (57 of 1,031 local eQTLs for BYxRM, 
#3 of 721 local eQTLs for YJM981xCBS2888, and 56 of 1,193 local eQTLs for YPS163xYJM145) had cell-cycle stage local eQTL interactions
#at a FDR of 5% across our mapping panels. We expanded our search to trans eQTLs and found that 790 (24.4%) of our 3,238 trans eQTLs
#had cell-cycle stage interactions (401 of 1,562 trans eQTLs for BYxRM, 195 of 1,126 trans eQTLs  for YJM981xCBS2888, 194 of 550 trans eQTLs
#for YPS163xYJM145). 
combined_objects$cis_table %>% summarise(has_int=sum(FDR < 0.05),n=n())
combined_objects$cis_table %>% group_by(cross) %>% summarise(has_int=sum(FDR < 0.05),n=n())
combined_objects$cis_table %>% summarise(has_int=sum(has_interaction),n=n())
combined_objects$cis_table %>% group_by(cross) %>% summarise(has_int=sum(has_interaction),n=n())
combined_objects$hotspot_peaks %>%  summarise(has_int=sum(has_interaction_trans),n=n())
combined_objects$hotspot_peaks %>% group_by(cross) %>%  summarise(has_int=sum(has_interaction_trans),n=n())


combined_objects$cis_table %>% filter(FDR < 0.05) %>% summarise(has_int=sum(has_interaction),n=n()-sum(has_interaction))

combined_objects$hotspot_peaks %>%  summarise(has_int=sum(has_interaction_trans),n=n()-sum(has_interaction_trans))
### ### ###Trans eQTLs were more likely than local eQTLs to have a cell-cycle stage interaction (OR=7.8, Fisher’s exact test, P < 2.2e-16), 
#which suggests that trans eQTLs depend more on the state of the cells than local eQTLs. 
mat1 = matrix(c(116,2829,790,2448),byrow = T,ncol=2,nrow=2)
mat1_res= fisher.test(mat1)
1/mat1_res$estimate


### We mapped a total of 20 unique cell-cycle occupancy QTLs (4 cell-cycle occupancy QTL for BYxRM, 6 cell-cycle occupancy for YJM981xCBS2888,
#and 10 cell-cycle occupancy QTL for YPS163xYJM145)

#combined_objects
#combined_objects$noise$NOISE_CIS_M

combined_objects$cell_cycle_lods %>% group_by(experiment) %>% summarise(n=n(),n_bins=length(unique(bin)))

table(cross_data$A$trans$cell_cycle_lods$bin)

table(cross_data$`3004`$trans$cell_cycle_lods$bin)

table(cross_data$B$trans$cell_cycle_lods$bin)

combined_objects$cis_eqtl_ase %>% filter(p_adj < 0.05) %>% group_by(cross) %>% summarise(has_int=sum(has_interaction),n=n())
