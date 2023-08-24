gpa_var = cross_data$B$trans$cell_cycle_lods[5,]$marker


p1 = cross_data$B$trans$cell_cycle_beta_all %>% filter(marker == gpa_var) %>% ggplot(aes(y=beta,x=cell_cycle)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=beta - 1.96*se,ymax=beta + 1.96*se),width=.2) + theme_classic() +
  ylab("Cell-cycle occupancy effect") + xlab("Cell cycle stage") + theme(text=element_text(size=18))# + ggtitle("GPA1")



gt = cross_data$B$trans$segdata$Gsub[,gpa_var]
af = cross_data$B$trans$segdata$cell.cycle.df
af$gt = gt
af$is_g1 = ifelse(af$cell_cycle=="G1",1,0)
per_cc_percent = af %>% group_by(is_g1) %>% summarise(af=mean(gt))



cross_data$B$trans$hotspot_enrichments_and_overlaps$`chrVIII:46887-140660_51`$enrich_mf %>% filter(type == "BP")%>% 
  ggplot(aes(x=(log2FC),y=reorder(Term,log2FC),fill=(classic)))+ 
  geom_bar(stat="identity") + theme_classic() + scale_fill_viridis_c() + theme(text=element_text(size=18)) #+ 
  #theme_
library(GO.db)


xx = as.list(cross_data$B$trans$hotspot_enrichments_and_overlaps$`chrVIII:46887-140660_51`$enrich_mf,num_char=1000)

GO = as.list(GOTERM)
terms = c()
for(term in xx$GO.ID){
terms=c(terms,GO[[term]]@Term)
  
}

cross_data$B$trans$hotspot_enrichments_and_overlaps$`chrVIII:46887-140660_51`$enrich_mf$Term2 = terms

p2 = cross_data$B$trans$hotspot_enrichments_and_overlaps$`chrVIII:46887-140660_51`$enrich_mf %>% filter(type == "BP")%>% 
  ggplot(aes(x=(log2FC),y=reorder(Term2,log2FC),fill=(classic)))+ 
  geom_bar(stat="identity") + theme_classic() + scale_fill_viridis_c() + theme(text=element_text(size=18)) + ylab("GO TERM") +
  xlab(expression(log[2]*" fold-change"))

(p1)  /(p2)                                                                                              
pg = cowplot::plot_grid(p1,p3,p2,p3,nrow=2,labels=c("A","C","B","D"),label_size = 18)

p3 = ggplot() + geom_point()

ggsave(pg,file="figure3.png")
