###

###
source("plot_fx.R")


sum(cross_data$`3004`$cis$cis_test_with_disp_expr$FDR < 0.05)
sum(cross_data$B$cis$cis_test_with_disp_expr$FDR < 0.05)

dim(cross_data$B$trans$hotspot_peaks)
dim(cross_data$`3004`$trans$hotspot_peaks)


length(table(cross_data$B$trans$hotspot_list$bin))
length(table(cross_data$`3004`$trans$hotspot_list$bin))

p3004_bar  = make_hotspot_hist(cross_data$`3004`$trans$hotspot_table,ylim=c(0,100))
pb_bar = make_hotspot_hist(cross_data$B$trans$hotspot_table,ylim=NULL)
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

pb = make_hotspot_map(cross_data$B$trans$combined_peaks,title = "YPS163 x YJM144")
p3004 = make_hotspot_map(cross_data$`3004`$trans$combined_peaks,title = "YJM981 x CBS2888")

#(p3004|pb)/(p3004_bar|pb_bar) + plot_annotation(tag_levels = "A")


(p3004 / p3004_bar) | (pb / pb_bar) + plot_annotation(tag_levels = "A")

#cowplot::plot_grid(pb,p3004,pb_bar,p3004_bar,nrow=2,labels=c("A","B","C","D"),label_size = 18,align="hv")
ggsave("figures/figure2.svg")
#p4 = ggplot(cross_data$B$trans$hotspot_table,aes(pos,count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
#  geom_bar(stat="identity", width=50000)+ geom_hline(yintercept = cross_data$B$trans$hotspot_threshold,color="red")  + 
#  xlab('Variant position (Mb)')+ylab('Distal eQTLs')+ scale_alpha(guide = 'none') + 
#  facet_grid(~chr,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
#  theme_bw()  + scale_x_continuous(breaks = c(.5e5,1e6,1.5e6,2e6),labels = function(x)x/1e6)


cross_data$A$trans$hotspot_list %>% group_by(bin) %>% summarise(pos=hotspot_pos[1])

#cross_data$
