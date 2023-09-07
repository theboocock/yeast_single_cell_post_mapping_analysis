
### Count local eQTLs
combined_objects$noise$ASE %>% filter(p_adj_ase < 0.05) %>% group_by(cross) %>% summarise(n=n())


combined_objects$noise$ASE_M


combined_objects$noise$ASE %>% filter(has_interaction) %>% group_by(cross) %>% summarise(n=n())
