p1 = combined_objects$noise$NOISE_WITH_REP %>% filter(p_adj_disp.ASE < 0.05| p_adj_ase.ASE < 0.05) %>% ggplot(aes(theta.ASE,theta.ASE_REP)) + geom_point() + 
  ylim(c(-5,5)) + xlim(c(-5,5)) + stat_cor(method="spearman") + theme_bw() + ylab(expression(Delta*ln*"(noise) discovery")) + 
  xlab(expression(Delta*ln*"(noise) replication"))
p2 = combined_objects$noise$NOISE_WITH_REP %>% filter(p_adj_disp.ASE < 0.05 | p_adj_ase.ASE < 0.05) %>% ggplot(aes(estimate.cond.ASE,estimate.cond.ASE_REP)) + geom_point() + 
  ylim(c(-5,5)) + xlim(c(-5,5)) + stat_cor(method="spearman") + theme_bw() + ylab(expression(Delta*ln*"(expression) discovery")) + 
  xlab(expression(Delta*ln*"(expression) replication"))
plot_grid(p1,p2)

#ab = combined_objects$noise$NOISE_WITH_REP %>% filter(p_adj_disp.ASE < 0.05)
