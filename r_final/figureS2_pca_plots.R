p_one = DimPlot(rds$cc_seurat, group.by="cell_cycle",label = T,label.size = 9,pt.size = 2,reduction="pca",dims = c(1,2))  + ggtitle("") + scale_color_brewer(palette = "Set2")
p_two = DimPlot(rds$cc_seurat, group.by="cell_cycle",label = T,label.size = 9,pt.size = 2,reduction="pca",dims = c(3,4))  + ggtitle("") + scale_color_brewer(palette = "Set2")
p_three = DimPlot(rds$cc_seurat, group.by="cell_cycle",label = T,label.size = 9,pt.size = 2,reduction="pca",dims = c(1,5)) + ggtitle("") + scale_color_brewer(palette = "Set2")
p_four= DimPlot(rds$cc_seurat, group.by="cell_cycle",label = T,label.size = 9,pt.size = 2,reduction="pca",dims = c(1,6)) + ggtitle("") + scale_color_brewer(palette = "Set2")
plot_grid(p_one, p_two, p_three, p_four, labels = c("A","B","C","D"))

ggsave("fig_final/s2_pca.png",width=16,height=12,dpi=300,device = png)
ggsave("fig_final//svg/s2_pca.svg",width=16,height=12,device = png)
