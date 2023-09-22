rds = readRDS("../rproj/out/cell_cycle/00_BYxRM_480MatA_1/cell_cycle_final.RDS")
rds2 = readRDS("../rproj/out/cell_cycle/00_BYxRM_480MatA_2//cell_cycle_final.RDS")


pcc = DimPlot(rds$cc_seurat,group.by ="cell_cycle",label = T,label.size = 18,pt.size = 2) + 
  xlab("UMAP 1") + ylab("UMAP 2") + theme(text=element_text(size=18)) + ggtitle("") 

pcc2=DimPlot(rds2$cc_seurat,group.by ="cell_cycle",label = T,label.size = 18,pt.size = 2) + 
  xlab("UMAP 1") + ylab("UMAP 2") + theme(text=element_text(size=18)) + ggtitle("")
pcc_hist= rds$cc_seurat@meta.data %>% group_by(cell_cycle) %>% summarise(n=n()) %>% summarise(cc=cell_cycle,prop=n/sum(n)) %>%
  ggplot(aes(y=prop,x=cc)) + geom_bar(stat="identity") + theme_bw()  + theme(text=element_text(size=18)) +
  ylab("Proportion of cells") + xlab("Cell-cycle stage") 
pcc_hist2= rds2$cc_seurat@meta.data %>% group_by(cell_cycle) %>% summarise(n=n()) %>% summarise(cc=cell_cycle,prop=n/sum(n)) %>%
  ggplot(aes(y=prop,x=cc)) + geom_bar(stat="identity") + theme_bw() + theme(text=element_text(size=18)) + 
  ylab("Proportion of cells") + xlab("Cell-cycle stage")+ ggtitle("Single-cell run 2")


anno$type_marker = paste(anno$Type,anno$Marker)

pf = FeaturePlot(rds$cc_seurat,anno$Marker,combine = F)  
pf = lapply(seq_along(pf), function(x){pf[[x]] + labs(title=anno$type_marker[x])})
#print(pf)
#plot_grid(p1,p2,pcc,pcc2)
pf1 = plot_grid(plotlist = pf,ncol=4)

pf = FeaturePlot(rds2$cc_seurat,anno$Marker,combine = F)  
pf = lapply(seq_along(pf), function(x){pf[[x]] + labs(title=anno$type_marker[x])})
#print(pf)
pf2 = plot_grid(plotlist = pf,ncol=4)


#df$cc
#t="identity")

plot_grid(pcc_hist +ggtitle("Single-cell run 1"),pcc_hist2,pcc,pcc2,labels=c("A","B","C","D"))
ggsave("fig_final//s2_cc_overview.png",bg="white",width=16,height=12,dpi=300)

