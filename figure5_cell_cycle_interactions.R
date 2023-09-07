### FIgure 3 cell-cycle interactionns #

#DimPlot(cross_data$A$cis$one_pot_seurat$`04_ByxRM_51-`$cc_in$cc_seurat, group.by="cell_cycle")
a_cis = cross_data$A$cis$cis_test_with_disp_expr %>% filter(FDR < 0.05)
table(ab$has_interaction)
cis = combined_objects$cis_table %>% filter(FDR < 0.05) %>% group_by(cross,has_interaction) %>% summarise(n=n()) %>% mutate(type="cis")
trans = combined_objects$hotspot_peaks %>% filter(FDR < 0.05) %>% group_by(cross,has_interaction_trans)  %>% rename(has_interaction=has_interaction_trans)%>% summarise(n=n()) %>% mutate(type="trans")
#rbind(trans,cis) %>% dplyr::count(type,has_interaction) %>% group_by(type) %>% summarise(frac=n/n())



p1 = rbind(trans,cis) %>% ggplot(aes(y=n,x=type,fill=has_interaction))  + geom_bar(stat="identity")
p2 = rbind(trans,cis) %>% ggplot(aes(y=n,x=type,fill=has_interaction))  + geom_bar(stat="identity") + facet_wrap(~cross)

combined_objects$hotspot_peaks %>% filter(FDR < 0.05) %>% group_by(has_interaction_trans) %>% summarise(n=n())
cross_data$A$trans$cell_cycle_beta_all$marker_idx = 1:nrow(cross_data$A$trans$cell_cycle_beta_all)


p1 = cross_data$A$trans$cell_cycle_beta_all %>% filter(chrom != "chrIII") %>%
  mutate(z=beta/se) %>% ggplot(aes(y=lod,x=pos,group=cell_cycle,color=cell_cycle))  + 
  geom_hline(yintercept = cross_data$B$trans$cell_cycle_threshold)+ geom_point() +geom_line(linewidth=1.5) + 
  facet_grid(~chrom_short_f,scales="free_x",space="free")+ylab("LOD") +
  scale_x_continuous(breaks = c(.5e6,1e6,1.5e6,2e6),labels = function(x)x/1e6) +
  #scale_y_continuous(breaks =c(.5e6,1e6,1.5e6,2e6), labels =  function(x)x/1e6) +
  xlab("Variant position (Mb)")  + ggtitle("BYxRM") + theme_classic() + scale_color_brewer(palette="Dark2")  + theme(text=element_text(size=18)) + ylim(c(0,18))

p2 = cross_data$`3004`$trans$cell_cycle_beta_all %>% filter(chrom != "chrIII") %>%
  mutate(z=beta/se) %>% ggplot(aes(y=lod,x=pos,group=cell_cycle,color=cell_cycle))  + 
  geom_hline(yintercept = cross_data$B$trans$cell_cycle_threshold)+ geom_point() +geom_line(linewidth=1.5) + 
  
  scale_x_continuous(breaks = c(.5e6,1e6,1.5e6,2e6),labels = function(x)x/1e6) +
  facet_grid(~chrom_short_f,scales="free_x",space="free") +ylab("LOD")  +xlab("Variant position (Mb)")+ ggtitle("YJM981xCBS2888") + theme_classic() + scale_color_brewer(palette="Dark2")  + theme(text=element_text(size=18))+ ylim(c(0,18))

p3 = cross_data$B$trans$cell_cycle_beta_all %>% filter(chrom != "chrIII") %>%
  mutate(z=beta/se) %>% ggplot(aes(y=lod,x=pos,group=cell_cycle,color=cell_cycle))  + 
  geom_hline(yintercept = cross_data$B$trans$cell_cycle_threshold)+ geom_point() +geom_line(linewidth=1.5) + 
  
  scale_x_continuous(breaks = c(.5e6,1e6,1.5e6,2e6),labels = function(x)x/1e6) +
  facet_grid(~chrom_short_f,scales="free_x",space="free") +ylab("LOD")  +xlab("Variant position (Mb)") + ggtitle("YPS163xYJM145") + theme_classic() + scale_color_brewer(palette="Dark2")  + theme(text=element_text(size=18)) + ylim(c(0,18))

#cross_data$B$trans$cell_cycle_beta_all %>% filter(chrom == "chrVIII") %>%
#  mutate(z=beta/se) %>% ggplot(aes(y=z,x=pos,group=cell_cycle,color=cell_cycle))  +  geom_point() +geom_line(linewidth=1.5) + 
#  facet_grid(~chrom,scales="free_x",space="free")  + ggtitle("B") + theme_classic() + scale_color_brewer(palette="Dark2")  + theme(text=element_text(size=18))


m1 = cross_data$A$trans$cell_cycle_lods[5,]$marker
p4 =cross_data$A$trans$cell_cycle_beta_all %>% filter(marker == m1) %>% ggplot(aes(y=beta,x=cell_cycle)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=beta - 1.96*se,ymax=beta + 1.96*se),width=.2) + theme_classic() +
  ylab("Cell-cycle occupancy effect") + xlab("Cell cycle stage") + theme(text=element_text(size=18)) + ggtitle(expression(italic("MKT1")))
  

for(i in 1:length(cross_data$B$trans$cell_cycle_lods)){
  #print()
  print(i)
  m1 =  cross_data$B$trans$cell_cycle_lods[i,]$marker
  
  print(cross_data$B$trans$cell_cycle_beta_all %>% filter(marker == m1) %>% ggplot(aes(y=beta,x=cell_cycle)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=beta - 1.96*se,ymax=beta + 1.96*se),width=.2) + theme_classic() +
    ylab("Cell-cycle occupancy effect") + xlab("Cell cycle stage") + theme(text=element_text(size=18)) + ggtitle(expression(italic("MKT1"))))
  
}


m3 = cross_data$`3004`$trans$cell_cycle_lods[3,]$marker
p5= cross_data$`3004`$trans$cell_cycle_beta_all %>% filter(marker == m3) %>% ggplot(aes(y=-beta,x=cell_cycle)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=-beta - 1.96*se,ymax=-beta + 1.96*se),width=.2) + theme_classic() +
  ylab("Cell-cycle occupancy effect") + xlab("Cell cycle stage") + theme(text=element_text(size=18)) + ggtitle(expression(italic("CYR1")))
m_gpa = cross_data$B$trans$cell_cycle_lods[5,]$marker
p6 = cross_data$B$trans$cell_cycle_beta_all %>% filter(marker == m_gpa) %>% ggplot(aes(y=-beta,x=cell_cycle)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=-beta - 1.96*se,ymax=-beta + 1.96*se),width=.2) + theme_classic() +
  ylab("Cell-cycle occupancy effect") + xlab("Cell cycle stage") + theme(text=element_text(size=18)) + ggtitle(expression(italic("GPA1")))


#p1 + p2 + p3 + p4 + p5
library(patchwork)
attempt_at_plot = ((p1/p2/p3) + plot_layout(guides="collect") | (p4/p5/p6)) + plot_layout(widths=c(2,.8))
attempt_at_plot
#m3= 
ggsave("figures/figure5_cell_cycle.svg",width=16,height=12)



#cell_cycle_beta_df 

gpa_var = cross_data$B$trans$cell_cycle_lods[5,]$marker

