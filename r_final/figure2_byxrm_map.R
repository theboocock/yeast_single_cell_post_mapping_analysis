source("plot_fx.R")
pcompare = ap_cis$hmmGenoCis %>% mutate(single_cell_sig = p_adj < 0.05) %>% 
  inner_join(ap_cis$prevGenoCis,by=c("gene")) %>% ggplot(aes(y=beta.y,x=beta.x,color=single_cell_sig)) + 
  geom_point(size=1.5) +
  geom_abline() + theme_bw()  + stat_cor(method="spearman",size=8,show.legend = F) +
  theme(text=element_text(size=18))  + xlab("local eQTL effect (single-cell HMM)") + 
  ylab("local eQTL effect (matched whole genome)") + #+ theme_ +
  scale_color_brewer(name="HMM FDR",labels=c(">= 0.05","< 0.05"),palette ="Dark2") +   
  guides(color = guide_legend(override.aes = list(size = 8)))
library(patchwork)
#p2 + (pcompare/p3) + plot_annotation(tag_levels = "A")


#p1=  magick::image_read_svg("../figure_and_tables_paper/Figures/one_pot.svg",width = 2000)
#p4 = ggdraw() + draw_image(p1)

#plot_grid(p4,nrow=2)

#p4 = ggdraw() + draw_image("../figure_and_tables_paper/Figures/one_pot.png")




by_new = cross_data$A$trans$hotspot_table
by_old = cross_data$A_bulk$trans$hotspot_table
#p_new  = make_hotspot_hist(by_new)
#p_new

p2 = ggplot(by_new,aes(pos,count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
  #geom_bar(stat="identity", width=50000)+
  geom_point(data=chrom_lengths_plot,aes(x=pos,y=0),alpha=0) +#
  xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
  facet_grid(~chrom_f,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme_classic() #+ # theme(panel.spacing=unit(0, "lines"))

pf= ggplot(by_new,aes(x=pos,y=count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
  geom_bar(stat="identity",width=50000)+
  geom_point(data=chrom_lengths_plot,aes(x=pos,y=0),alpha=0) +#
  xlab('')+ylab('')+ scale_alpha(guide = 'none') +
  facet_grid(~chrom_f,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme_classic() #+ # theme(panel.spacing=unit(0, 

chrom_old = "NA"
#pplot_grid(p2,pf,nrow=2)
for(i in 1:nrow(by_new)){
  pos = by_new$pos[i]
  chrom = by_new$chr[i]
  if(i==0){
    start = 0
  }else{
    start =by_new$end[i-1] + 1 #- 25000
    if(chrom_old == chrom){
      start = start
    }else{
      start = 0
    }
  }
  length = chrom_lengths$lengths[chrom_lengths$chrom == chrom]
  if(i==nrow(by_new)){
    end = length 
  }else{
    chrom_next = by_new$chr[i+1]
    end = by_new$pos[i+1] - 25000
    if(chrom_next == chrom){
      end = end
    }else{
      end = length
    }
  }
  #start = max(0,pos-25000)
  #end = min(length,pos+25000)  
  by_new$start[i] = start
  by_new$end[i] = end
  chrom_old = by_new$chr[i]
}


#apply(by_new,1,function(x){max(as.numeric(x[2]) - 25000,})

by_new = clean_up_hotspots_for_raster(by_new)
by_old = clean_up_hotspots_for_raster(by_old)

p2  = by_new %>% ggplot(aes(ymax=count, ymin=0,xmin=start,xmax=end)) + geom_rect() + facet_grid(~chrom_f,scales="free",space="free") + 
xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
  facet_grid(~chrom_f,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme_classic() #+  #theme(panel.spacing=unit
#+


p3  = by_old %>% ggplot(aes(ymax=count, ymin=0,xmin=start,xmax=end)) + geom_rect() + facet_grid(~chrom_f,scales="free",space="free") + 
  xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
  facet_grid(~chrom_f,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme_classic() #+  #theme(pa


#p3 = ggplot(by_old,aes(pos,count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
#  geom_bar(stat="identity", width=50000)+  geom_point(data=chrom_lengths_plot,aes(x=pos,y=0),alpha=0)+
#  xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
#  facet_grid(~chrom_f,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
#  theme_classic() #+  #theme(panel.spacing=unit(0, "lines"))

plot2 = p3 + geom_hline(yintercept = cross_data$A_bulk$trans$hotspot_threshold,color="red") +  theme(text=element_text(size=18),strip.background = element_blank(),
                                                                                                     strip.text.x = element_blank()) +#+  theme(plot.margin = unit(c(0,0,0,0), "cm"))  +
  scale_y_reverse() + xlab("Variant position (Mb)") + ylab("") + scale_x_continuous(breaks = c(.5e6,1e6,1.5e6,2e6),labels = function(x)x/1e6) + theme(panel.border = element_rect(fill = NA, 
                                                                                                                                                                                                colour = "grey20")) 
plot1 = p2 + geom_hline(yintercept = cross_data$A$trans$hotspot_threshold,color="red")  + theme(text=element_text(size=18),axis.title.x=element_blank(),
                                                                                                axis.text.x=element_blank(),
                                                                                                axis.ticks.x=element_blank(),strip.text.x = element_blank()) + theme(panel.border = element_rect(fill = NA, 
                                                                                                                                                                                                 colour = "grey20"))  
  #coord_cartesian(ylim=c(0,100))#+ scale_y_continuous(expand=c(0,10))
chrom_lengths_plot = chrom_lengths_plot %>% mutate(chrom_f = chrom_short_f)
#plot1 = plot1 + 
#plot2 = plot2 + geom_point(data=chrom_lengths_plot,aes(x=pos,y=0),alpha=0) + geom_bar(stat="identity")


pb = plot_grid(plot1,plot2,nrow=2,align = "v")

pa = make_hotspot_map(cross_data$A$trans$combined_peaks,title = "") + xlab("")
plot_grid(pa,plot1,plot2,align="v",axis="tblr",nrow=3,rel_heights = c(2,.4,.4))

#plot_grid(pa,(plot1/plot2),rel_heights = c(2,1),nrow=2,align="hv")
#(pa)/((plot1 )/(plot2)) + plot_layout(heights=c(2,.5))

ggsave("i/figure2_bulk_maps.svg",width=16,height=12)
ggsave("figures/figure2_bulk_maps.png",width=16,height=12)

#p4 / (pcompare+p3) + plot_annotation(tag_levels = "A") + plot_layout(height=c(1,2))
#ggsave("figures/figure1.png")

#pg = plot_grid(pcompare,p3,labels=c("B","C"),label_size = 18)
#plot_grid(p4,pg,labels=c("A",""),label_size = 18,nrow=2,rel_heights = c(1,2))
#ggsave("figures/figure1.png",dpi=300,bg = "white",width=8,height=6)

