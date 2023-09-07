combined_objects$noise$ASE  = combined_objects$noise$ASE %>% mutate(sig_new_filt =  (  combined_objects$noise$ASE$p_adj_ase > 0.05 | combined_objects$noise$ASE$estimate.cond*combined_objects$noise$ASE$estimate.disp < 0) & p_adj_disp < 0.05)# & p_adj_disp < 1e-3)# %>%
 # group_by(cross) %>% summarise(n=n())

sum(combined_objects$noise$ASE$sig_new_filt,na.rm=T)
combined_objects$noise$ASE  %>% group_by(cross,sig_new_filt) %>% summarise(n=n())
#combined_objects$noise$ASE  = combined_objects$noise$ASE %>% mutate(sig_new_filt =  (combined_objects$noise$ASE$estimate.cond*combined_objects$noise$ASE$estimate.disp < 0) & p_adj_disp < 0.05)# %>%
#combined_objects$noise$ASE %>% filter()
sum(combined_objects$noise$ASE$sig_new_filt,na.rm=T)
combined_objects$noise$ASE_EMMEANS=  combined_objects$noise$ASE_EMMEANS %>% mutate(theta=log(1/exp(emmean.disp))) 
big_drop = combined_objects$noise$ASE %>% filter(sig_new_filt) %>% filter(p_adj_disp < 1e-4) %>% filter(estimate.cond*estimate.disp< 0)# %>% ggplot(aes(x=abs(estimate.disp))) + geom_histogram()
big_drop = big_drop %>% group_by(cross) %>% slice_max(p_adj_disp,n=10)


#combined_objects$noise$ASE %>% filter(sig_new_filt) %>% filter(p_adj_disp < 1e-3) %>% ggplot(aes(x=abs(estimate.disp))) + geom_histogram()
big_drop_sub = big_drop %>% dplyr::select(gene,cross)

#big_drop %>% 


combined_objects$noise$ASE %>% filter(gene == "YGR192C")

segment_subset = combined_objects$noise$ASE_EMMEANS %>% pivot_wider(id_cols = c("gene","cross","cross2"),names_from = "geno",values_from = c("emmean.mean","theta","SE.mean","SE.disp")) %>%
  inner_join(big_drop_sub,by=c("gene","cross"))

segment_subset = segment_subset %>% inner_join(genes_name_trans,by=c("gene"="gene_id"))

segment_subset = segment_subset %>% mutate(avg=(emmean.mean_A + emmean.mean_B)/2,disp=(theta_A + theta_B)/2)

#segment_subset
combined_objects$noise$ASE %>% filter(sig_new_filt) 
segment_subset_long = segment_subset %>% pivot_longer(!gene & !cross & !gene_name) %>% filter(!grepl("SE",name))
#segment_subset_long %>% ifelse(grepl())


avg = segment_subset %>% dplyr::select(gene,cross,gene_name,avg,disp) %>% mutate(allele=NA)
not_avgA = segment_subset %>% dplyr::select(gene,cross,gene_name, emmean.mean_A,theta_A) %>% mutate(allele="A")
not_avgB = segment_subset %>% dplyr::select(gene, cross,gene_name,emmean.mean_B,theta_B) %>% mutate(allele="B")

col_names_str = c("gene","cross","gene_name","mean","disp","allele")
colnames(avg) = col_names_str
colnames(not_avgA) = col_names_str
colnames(not_avgB) = col_names_str

avg_m = rbind(avg %>% mutate(label=gene_name),not_avgA %>% mutate(label=""),not_avgB %>% mutate(label=""))# %>% mutate(ifelse())
library(ggrepel)

avg_m = avg_m %>% mutate(cross2  = case_when(
  cross == "3004" ~ "YJM981 (red) x CBS2888 (purple)",
  cross == "B" ~ "YPS163 (red) x YJM145 (purple)",
  cross == "A" ~ "BY (red) x RM (purple)"
))
#segment_subset %>% 
#combined_objects$noise$ASE_EMMEANS %>%  ggplot(aes(y=theta,x=emmean.mean)) + geom_point(alpha=0.1) + geom_segment(data=segment_subset,aes(x=emmean.mean_A,y=theta_A,xend=emmean.mean_B,yend=theta_B),color="red",size=.5)+ 
# geom_label_repel(data=avg_m,aes(x=mean,y=disp,label=label),box.padding = 2,max.overlaps = Inf) +
#  geom_point(data=segment_subset,aes(x=emmean.mean_A,y=theta_A),color="red") + geom_point(data=segment_subset,aes(x=emmean.mean_B,y=theta_B),color="red") + 
#  geom_errorbar(data=segment_subset,aes(x=emmean.mean_A,y=theta_A,xmin=emmean.mean_A - 1.96*SE.mean_A,xmax=emmean.mean_A + 1.96*SE.mean_A,
#                                        ymin=theta_A - 1.96*SE.disp_A,ymax=theta_A + 1.96*SE.disp_A),color="red") + 
#  geom_errorbarh(data=segment_subset,aes(x=emmean.mean_A,y=theta_A,xmin=emmean.mean_A - 1.96*SE.mean_A,xmax=emmean.mean_A + 1.96*SE.mean_A,
#                                        ymin=theta_A - 1.96*SE.disp_A,ymax=theta_A + 1.96*SE.disp_A),color="red") + 
#  geom_point(data=segment_subset,aes(x=emmean.mean_B,y=theta_B),color="red") + geom_point(data=segment_subset,aes(x=emmean.mean_B,y=theta_B),color="red") + 
#  geom_errorbar(data=segment_subset,aes(x=emmean.mean_B,y=theta_B,xmin=emmean.mean_B - 1.96*SE.mean_A,xmax=emmean.mean_B + 1.96*SE.mean_B,#
#                                        ymin=theta_B - 1.96*SE.disp_B,ymax=theta_B + 1.96*SE.disp_B),color="red") +
#  geom_errorbarh(data=segment_subset,aes(x=emmean.mean_B,y=theta_B,xmin=emmean.mean_B - 1.96*SE.mean_A,xmax=emmean.mean_B + 1.96*SE.mean_B,
#                                        ymin=theta_B - 1.96*SE.disp_B,ymax=theta_B + 1.96*SE.disp_B),color="red") +
#  facet_wrap(~cross2)  + stat_smooth(level=.99) + theme_classic()  + xlim(c(-7,5)) + ylim((c(-5,5))) + theme(text=element_text(size=18)) + xlab(expression(ln*"(expression)")) + ylab(expression(ln(1/theta)))
combined_objects$noise$ASE_EMMEANS %>%  ggplot(aes(y=theta,x=emmean.mean)) + geom_point(alpha=0.1) + geom_segment(data=segment_subset,aes(x=emmean.mean_A,y=theta_A,xend=emmean.mean_B,yend=theta_B),size=.5) + stat_smooth(level=.99) +
  geom_point(data=segment_subset,aes(x=emmean.mean_A,y=theta_A),color="red",size=5) +  #+ geom_point(data=segment_subset,aes(x=emmean.mean_B,y=theta_B),color="purple",size=2) + 
  geom_point(data=segment_subset,aes(x=emmean.mean_B,y=theta_B),color="purple",size=5) +
  geom_label_repel(data=avg_m,aes(x=mean,y=disp,label=label),box.padding = 2,max.overlaps = Inf,size=6) +# + geom_point(data=segment_subset,aes(x=emmean.mean_B,y=theta_B),color="purple",size=2) + 
  facet_wrap(~cross2) + theme_classic()  + theme(text=element_text(size=18)) + xlab(expression(ln*"(expression)")) + ylab(expression(ln(dispersion))) + coord_cartesian(ylim=c(-3,5),xlim=c(-5.5,5))



combined_objects$noise$ASE  %>% filter(gene == "YMR321C")
segment_subset  %>% filter(gene == "YMR321C")

ggsave(filename = "figures/figure4_noise.png",width=16,height=12)
ggsave(filename = "figures/figure4_noise.svg",width=16,height=12)

#segment_subset %>% filter(gene_name == "HSP12") %>% ggplot()

#segment_subset %>% 


## Haploid direction
## diploid direction





combined_objects$noise$ASE_EMMEANS  %>% filter(gene == "YFL014W") %>%ggplot(aes(y=theta,x=emmean.mean, color=geno)) + geom_point() + facet_wrap(~cross) 

#%>% ggplot(aes(y=))

#plot(combined_objects$noise$ASE$estimate.cond,combined_objects$noise$ASE$estimate)




#combined_objects$
  
  
rep_ase = combined_objects$noise$ASE_M=combined_objects$noise$ASE_M %>% filter(sig_trend_filter)
sum(rep_ase$estimate.disp.ASE * rep_ase$estimate.disp.REP > 0,na.rm=T)


sum(rep_ase$estimate.disp.ASE * rep_ase$estimate.disp.REP > 0 & rep_ase$p.value.disp.REP < 0.05,na.rm=T)

370/(length(na.omit(rep_ase$p.value.disp.REP)))


combined_objects$noise$ASE

combined_objects$noise$ASE %>% filter(p_adj_ase > .05 | (combined_objects$noise$ASE$estimate.cond*combined_objects$noise$ASE$estimate.disp < 0)) %>% filter(p_adj_disp < 0.05) %>% arrange(p_adj_disp
                                                                                                                                                                                           )


#combined_objects$noise$ASE$estimate.disp combined_objects$noise$ASE$estimate.disp


#combined_objects$noise$NOISE_CIS_M %>% ggplot(aes(x=estimate,y=Beta)) + geom_point() + facet_wrap(~cross.ASE)



combined_objects$cis_eqtl_ase_noise %>% ggplot(aes(x=Beta,y=estimate,color=FDR < 0.05)) + geom_point() + facet_wrap(~cross2,nrow=3) + theme_bw() + theme(text=element_text(size=18)) +
  coord_cartesian(ylim=c(-2,2)) +  stat_cor(method="spearman",label.x = c(-2),label.y=c(2,1.8),size=8)   + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + xlab("local eQTL effect (one-pot)")  + ylab("allele-specific expression effect") + 
  scale_color_brewer(name="One-pot FDR",labels=c(">= 0.05","< 0.05"),palette ="Dark2") +   
  guides(color = guide_legend(override.aes = list(size = 8)))#+ geom_hline() +




#whit= read.csv("~/whitkopp.csv")
#whit[grep("A|NN",whit$STRAIN),]
#whit2= read.csv("~/whitkopp2.csv")
#whit2 = whit2 %>% dplyr::rename(STRAIN = X.STRAIN.....A....STRAIN.....Z....STRAIN.....WT.Y1.1....STRAIN.....ZZ....STRAIN.....AA.)
#p1 = whit2 %>% filter(STRAIN == "A" | STRAIN == "Z") %>% ggplot(aes(y=YFP.MEAN.EFFECT,x=STRAIN)) + geom_boxplot()  + geom_jitter()
#p2 = whit %>% filter(STRAIN == "A" | STRAIN == "Z" ) %>% ggplot(aes(y=YFP.MEAN.EFFECT,x=STRAIN)) + geom_boxplot()  + geom_jitter()
#cowplot::plot_grid(p1,p2,nrow=2)
#whit2%>% filter(STRAIN == "A" | STRAIN == "Z" | STRAIN == "WT.Y1.1" | STRAIN == "ZZ" | STRAIN == "AA") %>% ggplot(aes(y=YFP.SD.CORRECT,x=STRAIN)) + geom_boxplot()
#whit %>% ggplot(aes(y=YFP.SD.CORRECT,x=STRAIN)) + geom_boxplot()
