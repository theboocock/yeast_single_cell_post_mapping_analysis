pcompare = ap_cis$hmmGenoCis %>% mutate(single_cell_sig = p_adj < 0.05) %>% 
  inner_join(ap_cis$prevGenoCis,by=c("gene")) %>% ggplot(aes(y=beta.y,x=beta.x)) + 
  geom_point(size=1.5) +
  geom_abline() + theme_bw()  + stat_cor(method="spearman",size=8,show.legend = F,cor.coef.name = "rho") +
  theme(text=element_text(size=24))  + xlab("local eQTL effect (HMM)") + 
  ylab("local eQTL effect (whole genome sequencing)") + #+ theme_ +
  scale_color_brewer(name="HMM FDR",labels=c(">= 0.05","< 0.05"),palette ="Dark2") + # + 
  guides(color = guide_legend(override.aes = list(size = 8)))
pcompare
ggsave("fig_final/s8.png",dpi=300,height = 12,width=16)
ggsave("fig_final/svg//s8.svg",height = 12,width=16)
