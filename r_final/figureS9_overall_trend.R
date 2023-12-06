ab = combined_objects$noise$ASE %>% mutate(is_sig = p_adj_disp < 0.05 | p_adj_ase < 0.05)  %>% mutate(not_sig = p_adj_disp > 0.05 & p_adj_ase > 0.05)

p2 = ab %>% filter(!is.na(is_sig)) %>% ggplot(aes(y=theta,x=estimate.cond,color=is_sig)) + geom_point() + ylim(-5,5) + xlim(-5,5)   +
  theme_bw()  + ylab(expression(Delta * ln(noise))) +
  xlab(expression(Delta * ln(expression))) + theme(text=element_text(size=18)) + scale_color_manual(name="Significant (FDR < 0.05)\nallele-specific effect",values = c("grey70","#4daf4a"),
                                                                                                    labels=c("No","Yes"))

cor(ab$estimate.cond,ab$estimate.disp,method="spearman",use="pairwise.complete.obs")

cor(ab$estimate.cond[!ab$is_sig],ab$estimate.disp[!ab$is_sig],method="spearman",use="pairwise.complete.obs")
cor(ab$estimate.cond[ab$is_sig],ab$estimate.disp[ab$is_sig],method="spearman",use="pairwise.complete.obs")

ab  = ab %>% mutate(gene_cross=paste(gene,"_",cross,sep=""))

combined_objects$noise$ASE_EMMEANS = combined_objects$noise$ASE_EMMEANS %>% mutate(gene_cross=paste(gene,"_",cross,sep=""))

em_subset = combined_objects$noise$ASE_EMMEANS[(combined_objects$noise$ASE_EMMEANS$gene_cross %in% ab$gene_cross[!is.na(ab$is_sig)]),]

p3 = em_subset %>%  ggplot(aes(y=theta,x=emmean.mean)) + geom_point(color="grey70") + xlim(-7,5) + ylim(-5,7) + theme_bw() + ylab(expression(ln(noise))) + xlab(expression(ln(expression))) +
  theme(text=element_text(size=18))

cor(em_subset$theta,em_subset$emmean.mean,method="spearman")

plot_grid(p3,p2,labels =c("A","B"), label_size = 18,rel_widths = c(.7,1))#,rel_widths = c(0.7,1))


ggsave("fig_final/svg/s9.svg",width=16,height=12)
ggsave("fig_final/s9.png",dpi=300,width=16,height=12,device=png)
# 142 rows removed out of 
# + stat_cor(method="spearman")                                   
