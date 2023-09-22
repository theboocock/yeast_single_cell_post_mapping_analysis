#We found a total of 20 unique cell-cycle occupancy QTLs in the three crosses
#(4 for cross A, 10 for cross B, and 6 for cross C; Figure 2A and 5A, Table S11)

combined_objects$cell_cycle_lods %>% group_by(experiment) %>% summarise(n=n(),n_bins=length(unique(bin)))
