
big_drop = combined_objects$noise$ASE %>% filter(sig_new_filt)# %>% filter(p_adj_disp < 1e-4) %>% filter(estimate.cond*estimate.disp< 0)# %>% ggplot(aes(x=abs(estimate.disp))) + geom_histogram()

sum(combined_objects$noise$ASE$p_adj_disp < 0.05,na.rm=T)

big_drop = big_drop %>% group_by(cross) %>% slice_min(p_adj_disp,n=10,with_ties = F)


big_drop %>% filter(gene == "YFL014W")
big_drop_sub = big_drop %>% dplyr::select(gene,cross)
segment_subset = combined_objects$noise$ASE_EMMEANS %>% pivot_wider(id_cols = c("gene","cross","cross2"),names_from = "geno",values_from = c("emmean.mean","theta","SE.mean","SE.disp")) %>%
  inner_join(big_drop_sub,by=c("gene","cross"))
segment_subset = segment_subset %>% inner_join(genes_name_trans,by=c("gene"="gene_id"))
segment_subset = segment_subset %>% mutate(avg=(emmean.mean_A + emmean.mean_B)/2,disp=(theta_A + theta_B)/2)
avg = segment_subset %>% dplyr::select(gene,cross,gene_name,avg,disp) %>% mutate(allele=NA)
not_avgA = segment_subset %>% dplyr::select(gene,cross,gene_name, emmean.mean_A,theta_A) %>% mutate(allele="A")
not_avgB = segment_subset %>% dplyr::select(gene, cross,gene_name,emmean.mean_B,theta_B) %>% mutate(allele="B")
col_names_str = c("gene","cross","gene_name","mean","disp","allele")
colnames(avg) = col_names_str
colnames(not_avgA) = col_names_str
colnames(not_avgB) = col_names_str
avg_m = rbind(avg %>% mutate(label=gene_name),not_avgA %>% mutate(label=""),not_avgB %>% mutate(label=""))# %>% mutate(ifelse())

avg_m = avg_m %>% mutate(cross2  = case_when(
  cross == "3004" ~ "CBS2888 (red) x YJM981 (purple)",
  cross == "B" ~ "YJM145 (red) x YPS163 (purple)",
  cross == "A" ~ "BY (red) x RM (purple)"
))
library(ggrepel)
avg_m %>% filter(gene_name == "HSP12")


nrow(combined_objects$noise$ASE_EMMEANS %>% dplyr::filter(theta > -3 & theta < 5) %>% dplyr::filter(emmean.mean > -5.5 & emmean.mean < 5))

nrow(combined_objects$noise$ASE_EMMEANS)

combined_objects$noise$ASE_EMMEANS %>%  ggplot(aes(y=theta,x=emmean.mean)) + geom_point(alpha=0.05) + geom_segment(data=segment_subset,aes(x=emmean.mean_A,y=theta_A,xend=emmean.mean_B,yend=theta_B),size=.5) + stat_smooth(level=.99) +
  geom_point(data=segment_subset,aes(x=emmean.mean_A,y=theta_A),color="red",size=5) +  #+ geom_point(data=segment_subset,aes(x=emmean.mean_B,y=theta_B),color="purple",size=2) + 
  geom_point(data=segment_subset,aes(x=emmean.mean_B,y=theta_B),color="purple",size=5) +
  geom_label_repel(data=avg_m,aes(x=mean,y=disp,label=label),box.padding = 2,max.overlaps = Inf,size=6) +# + geom_point(data=segment_subset,aes(x=emmean.mean_B,y=theta_B),color="purple",size=2) + 
  facet_wrap(~cross2) + theme_classic()  + theme(text=element_text(size=18)) + xlab(expression(ln*"(expression)")) + ylab(expression(ln(dispersion))) + coord_cartesian(ylim=c(-3,5),xlim=c(-5.5,5))
#ggsave("figure_fin  3.png",width=16,height=12)
ggsave("fig_final//4.svg",width=16,height=12)


cor(combined_objects$noise$ASE_EMMEANS$theta,combined_objects$noise$ASE_EMMEANS$emmean.mean,method="spearman")
