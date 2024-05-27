rds$cc_seurat@meta.data$cell_cycle = factor(rds$cc_seurat@meta.data$cell_cycle,levels=c("M/G1","G1","G1/S","S","G2/M"))
DoHeatmap(rds$cc_seurat[anno$Marker,], features = anno$Marker, group.by="cell_cycle",label = T)  + ggtitle("") + 
  scale_color_brewer(palette = "Set2")
ggsave("fig_final/s3_heatmap.png",width=16,height=12,dpi=300,device = png)
ggsave("fig_final//svg/s3_heatmap.svg",width=16,height=12,device = png)




