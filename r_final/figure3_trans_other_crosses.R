###

###




p3004_bar  = make_hotspot_hist(cross_data$`3004`$trans$hotspot_table,ylim=c(0,100))# + theme(strip.text.x = element_blank()) + theme(pane)
pb_bar = make_hotspot_hist(cross_data$B$trans$hotspot_table,ylim=NULL)# + theme(strip.text.x = element_blank())
#p3004_bar= ggplot(cross_data$`3004`$trans$hotspot_table,aes(pos,count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
#  geom_bar(stat="identity", width=50000)+geom_hline(yintercept = cross_data$`3004`$trans$hotspot_threshold,color="red")  + 
#  xlab('Variant position (Mb)')+ylab("Distal eQTLs")+ scale_alpha(guide = 'none') + 
#  facet_grid(~chrom_f,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
#  theme_classic() + scale_x_continuous(breaks = c(.5e6,1e6,1.5e6,2e6),labels = function(x)x/1e6) + coord_cartesian(ylim=(c(0,100))) +
#  theme(text=element_text(size=18))
#pb_bar= ggplot(cross_data$B$trans$hotspot_table,aes(pos,count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
#  geom_bar(stat="identity", width=50000)+geom_hline(yintercept = cross_data$B$trans$hotspot_threshold,color="red")  + 
#  xlab('Variant position (Mb)')+ylab("Distal eQTLs")+ scale_alpha(guide = 'none') + 
#  facet_grid(~chrom_f,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
#  theme_classic() + scale_x_continuous(breaks = c(.5e6,1e6,1.5e6,2e6),labels = function(x)x/1e6) +#
#  theme(text=element_text(size=18))

#geom_hline(yintercept = cross_data$A$trans$hotspot_threshold,color="red")   + theme_classic() +
#  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +  
##  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + scale_y_continuous(expand=c(0,0)) +
#  theme(text=element_text(size=18)) + coord_cartesian(ylim=c(0,100)) + ylab("Distal eQTLs (single-cell)

#+ scale_log10 +)

# + scale_x_continuous(breaks = c(5e6,1e7,1.5e7,2e7),labels = function(x)x/1e6)

library(patchwork)

#plot1 /plot2

pb = make_hotspot_map(cross_data$B$trans$combined_peaks,title = "YJM145 x YPS163 (Cross B)")
p3004 = make_hotspot_map(cross_data$`3004`$trans$combined_peaks,title = "CBS2888 x YJM981 (Cross C)")

plot_grid(pb,p3004,pb_bar,p3004_bar,nrow=2,labels = c("A","B","C","D"),label_size = 18,rel_heights = c(2,.6),align="hv",axis="lrtb")
#(p3004|pb)/(p3004_bar|pb_bar) + plot_annotation(tag_levels = "A")
#(p3004 / p3004_bar) | (pb / pb_bar) + plot_annotation(tag_levels = "A")

#cowplot::plot_grid(pb,p3004,pb_bar,p3004_bar,nrow=2,labels=c("A","B","C","D"),label_size = 18,align="hv")

ggsave("fig_final/main/staging/figure3.svg",width=26,height=16)#,width=26,height=16)

#2.6/2

#1.3 / 2
#30.1/()
#ggsave("figure_drafts//3.svg",width=26,height=26*.65)

#ggsave("figure_drafts//3.svg",width=26,height=16)
#p4 = ggplot(cross_data$B$trans$hotspot_table,aes(pos,count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
#  geom_bar(stat="identity", width=50000)+ geom_hline(yintercept = cross_data$B$trans$hotspot_threshold,color="red")  + 
#  xlab('Variant position (Mb)')+ylab('Distal eQTLs')+ scale_alpha(guide = 'none') + 
#  facet_grid(~chr,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
#  theme_bw()  + scale_x_continuous(breaks = c(.5e5,1e6,1.5e6,2e6),labels = function(x)x/1e6)


#cross_data$A$trans$hotspot_list %>% group_by(bin) %>% summarise(pos=hotspot_pos[1])

#cross_data$
