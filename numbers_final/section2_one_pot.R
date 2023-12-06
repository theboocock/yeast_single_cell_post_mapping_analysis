#We sampled ~250,000 MATa segregants and used scRNA-seq to
#measure the expression of 5,435 transcripts in 27,744 single cells (Table S4).
summary_table[2,]
#We mapped 1,031 local eQTLs at an FDR of 5%. We compared these results to those from bulk RNA-seq and genotyping in the same cross (Albert et al. 2018) and found that 717 (69.5%) of the 1,031 local eQTLs were also detected as statistically significant in that study, with an additional 108 local eQTLs showing effects in the same direction (Figures 1C and 2A; Table S7). Thus, 825 (80%) of
#the local eQTLs detected with the one-pot approach were supported by the bulk results

sum(one_pot_bulk$FDR < 0.05,na.rm=T)
717/1031
sum(one_pot_bulk$FDR < 0.05 &
      one_pot_bulk$Beta * one_pot_bulk$coefficient > 0 &one_pot_bulk$p_adj < 0.05,na.rm=T )

sum(one_pot_bulk$FDR < 0.05 & one_pot_bulk$Beta * one_pot_bulk$coefficient > 0,na.rm=T)
825/1031
#We mapped 1,562 distal eQTLs at an FDR of 5% (Figure 2A; Table S8).
nrow(cross_data$A$trans$hotspot_peaks_new)




A_hot
A_old_hot

A_hot %>% join_overlap_left(A_old_hot)
m_hotspots= A_hot %>% join_overlap_left(A_old_hot)
m_hotspots_old = A_old_hot %>% join_overlap_left(A_hot)
m_hotspots[!is.na(m_hotspots$bin.y),]
m_hotspots[is.na(m_hotspots$bin.y),]

m_hotspots_old[is.na(m_hotspots_old$hotspot_pos.y)]

# New hotspots


# We took advantage of the convenience of one-pot eQTL mapping and applied it to two additional yeast crosses, one between a clinical 
#strain (YJM145) and a soil strain (YPS163) both isolated in the United States

summary_table[3,]
joseph_annotation[grep("YJM145_b",joseph_annotation$Isolate.name),]
joseph_annotation[grep("YPS163_1b",joseph_annotation$Isolate.name),]

summary_table[4,]
joseph_annotation[grep("YJM981",joseph_annotation$Isolate.name),]
joseph_annotation[grep("YPS163_1b",joseph_annotation$Isolate.name),]

#We mapped a total of 1,914 local eQTLs in the new crosses (1,193 in cross B and 721 in cross C; Table S7), 
#as well as 1,626 distal eQTLs (550 in cross B and 1,126 in cross C; Figure 3A-B; Table S8)
nrow(cross_data$B$cis$cis_test_with_disp_expr %>% filter(FDR < 0.05) )
nrow(cross_data$`3004`$cis$cis_test_with_disp_expr %>% filter(FDR < 0.05))
#as well as 1,626 distal eQTLs (550 in cross B and 1,126 in cross C; Figure 3A-B; Table S8).
nrow(cross_data$B$trans$hotspot_peaks_new)
nrow(cross_data$`3004`$trans$hotspot_peaks_new)

# These distal eQTLs clustered into 13 hotspots (6 in cross C and 7 in cross B; Figure 3C-D). #
B_hot
C_hot

length(c(B_hot,C_hot,A_hot))

AC = A_hot %>% join_overlap_left(B_hot)
AC = AC[is.na(AC$hotspot_pos.y)] %>% join_overlap_left(C_hot)
AC[is.na(AC$hotspot_pos)]
# 7/12
BC = B_hot %>% join_overlap_left(C_hot)
BC = BC[is.na(BC$hotspot_pos.y)] %>% join_overlap_left(A_hot)
BC[is.na(BC$hotspot_pos)]
# 2/12
CA = C_hot %>% join_overlap_left(B_hot)
CA = CA[is.na(CA$hotspot_pos.y)] %>% join_overlap_left(A_hot)
CA[is.na(CA$hotspot_pos)]

