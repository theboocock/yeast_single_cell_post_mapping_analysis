source("pipeline/load_cc.R")
#source("pipeline/utils.R")
source("cell_cycle_annotation/cell_cycle_annotation_fx.R")
# Let's start with the 480 BYxRM #
library(tidyverse)
library(glue)
library(cowplot)
# Haploids
input_files = read.csv("data/single_cell_runs.csv")
input_files = input_files %>% filter(type == "haploids")
for(j in 1:nrow(input_files)){
  in_folder = input_files$name[j]
  print(in_folder)
  out_folder= glue("data/cell_cycle/out/{in_folder}")
  dir.create(out_folder)
  cell_filter = readRDS(glue("{out_folder}/{cellfilter_rds}"))
  out_png = glue("{out_folder}/feature_plot.png")
  png(out_png,width=1600,height = 800)
  p1 = FeaturePlot(cc$cc_seurat,anno$Marker[anno$Type == "M/G1"]) 
  p2 =  FeaturePlot(cc$cc_seurat,anno$Marker[anno$Type == "G1"])
  p3 = FeaturePlot(cc$cc_seurat,anno$Marker[anno$Type == "S"]) 
  p4 = FeaturePlot(cc$cc_seurat,anno$Marker[anno$Type == "G2/M"]) 
  print(plot_grid(p1,p2,p3,p4,labels = c("M/G1","G1","S","G2/M")))
  dev.off()
  #drive_upload(out_png, path=cell_cycle_folder,name="feature_plot.png")
  
  top10 <- cc$cc2 %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  out_png = glue("{out_folder}/top_genes_clusters.png")
  png(out_png,width=1600, height=800)
  print(DoHeatmap(cc$cc_seurat, features = top10$gene) + NoLegend())
  dev.off()
  out_png = glue("{out_folder}/clustering.png")
  png(out_png,width=1600,height = 800)
  print(DimPlot(cc$cc_seurat,label.size = 12,label=T,pt.size = 2))
  dev.off()
  #drive_upload(out_png, path=cell_cycle_folder,name="clustering.png")i
  cc = cellcycle_annotation(glue("{out_folder}/"),cell_filter)
  
  out_rds_cell_cycle = glue("{out_folder}/cell_cycle.RDS")
  saveRDS(cc,file=out_rds_cell_cycle)
  # out_all_rds_cell_cycle = glue("{out_folder}/all.RDS")
  # saveRDS(all,file=out_all_rds_cell_cycle)
  #drive_upload("tmp/cell_cycle.RDS",path=cell_cycle_folder,name="cell_cycle.RDS",type = "application/zip")
  all_data = read_and_normalize(glue("{out_folder}/"),cell_filter)
  ### ###
  
  out_rds_all =  glue("{out_folder}/all.RDS")
  saveRDS(all_data,file=out_rds_all)
  out_png = glue("{out_folder}/all_clustering.png")
  png(out_png,width=1600,height = 800)
  print(DimPlot(all_data$seurat_obj,label.size = 12,label=T,pt.size = 2))
  dev.off()
  
  top10 <- all_data$markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  out_png = glue("{out_folder}/all_top_genes_clusters.png")
  png(out_png,width=1600, height=800)
  print(DoHeatmap(all_data$seurat_obj, features = top10$gene) + NoLegend())
  dev.off()
  #drive_upload("tmp/all_seurat.RDS",path=cell_cycle_folder,name="all_seurat.RDS",type = "application/zip")
}
# Diploids
input_files = read.csv("data/single_cell_runs.csv")
input_files = input_files %>% filter(type == "diploids" & good == "Y")

for(j in 1:nrow(input_files)){
  in_folder = input_files$name[j]
  out_folder= glue("data/cell_cycle/out/{in_folder}")
  dir.create(out_folder)
  diploid_assignments = "diploid_assignments.RDS"
  diploid_assignments = readRDS(glue("{out_folder}/{diploid_assignments}"))
  in_seurat = glue("{out_folder}/")
  diploid_annotation(in_seurat = glue("{out_folder}/"),diploid_assignments = diploid_assignments, out_folder=out_folder)
}
