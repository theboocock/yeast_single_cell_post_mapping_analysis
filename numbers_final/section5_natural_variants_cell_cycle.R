#We found a total of 20 unique cell-cycle occupancy QTLs in the three crosses
#(4 for cross A, 10 for cross B, and 6 for cross C; Figure 2A and 5A, Table S11)

combined_objects$cell_cycle_lods %>% group_by(experiment) %>% summarise(n=n(),n_bins=length(unique(bin)))


combined_objects$cell_cycle_lod %>% filter(experiment == "A")


cross_data$`3004`$trans$hotspot_enrichments_and_overlaps$`chrX:397734-497167_224`$assoc_round_robin %>% filter(Beta < 0 & q.value < 0.05)
