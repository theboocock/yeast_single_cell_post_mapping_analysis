make_hotspot_hist = function(hotspot_table,ylim=NULL){
  p1= ggplot(hotspot_table,aes(pos,count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
    geom_bar(stat="identity", width=50000)+geom_hline(yintercept = cross_data$`3004`$trans$hotspot_threshold,color="red")  + 
    xlab('Variant position (Mb)')+ylab("Distal eQTLs")+ scale_alpha(guide = 'none') + 
    facet_grid(~chrom_f,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    theme_classic() + scale_x_continuous(breaks = c(.5e6,1e6,1.5e6,2e6),labels = function(x)x/1e6) + #coord_cartesian(ylim=(c(0,100))) +
    theme(text=element_text(size=18))
    if(!is.null(ylim)){
      print("HRE")
      p1 = p1 + coord_cartesian(ylim=ylim)
    }
  return(p1)
}

make_hotspot_map = function(combined_peaks,title="",fdr_filt = 0.05){
  p1  = combined_peaks %>% filter(FDR < fdr_filt)  %>% 
  sample_frac() %>% filter(tchr!= "chrmt") %>%
  ggplot(aes(x=pos,y=tpos)) +#+ geom_point(data=chrom_lengths_plot,aes(y=pos,x=pos))
  geom_point() +  geom_point(data=chrom_lengths_plot,aes(y=pos,x=pos),alpha=0) + 
  facet_grid(tchrom_short_f~chrom_short_f,scales="free",space = "free",switch="both") + 
  theme_classic() + scale_color_brewer(palette ="Dark2") + theme(text=element_text(size=18)) + 
  ylab("Transcript position (Mb)") + xlab("Variant position (Mb)") + 
  scale_x_continuous(breaks = c(.5e6,1e6,1.5e6,2e6),labels = function(x)x/1e6) +
  scale_y_continuous(breaks =c(.5e6,1e6,1.5e6,2e6), labels =  function(x)x/1e6) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "grey20")) 
  
  if(title != ""){
    p1 = p1 + ggtitle(title)
  }
 # + ggtitle("YPS163 x YJM454")
  return(p1)
}
