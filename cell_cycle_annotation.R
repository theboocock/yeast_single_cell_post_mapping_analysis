source("load_cc.R")
source("utils.R")
source("cell_cycle_annotation_fx.R")
# Let's start with the 480 BYxRM #
library(tidyverse)
library(glue)
library(googledrive)
library(glmnet)
library(cowplot)
# Note: Josh filters cells with > 10,000 UMIs. 
#in_mm = Sys.glob("/HOFFMAN/old_alignments/processed/processed/*/*outs/filtered_feature_bc_matrix/")
input_files = read.csv("single_cell_runs.csv")
input_files = input_files %>% filter(type == "haploids")
folder_in  = drive_find(q = "sharedWithMe and name = 'yeast'",n_max = 25)
results=drive_ls(folder_in)[drive_ls(folder_in)[,1] == "results",]
processed =drive_ls(folder_in)[drive_ls(folder_in)[,1] == "processed",]
results_files = drive_ls(results)
processed_files = drive_ls(processed)
#processed_files %>% select(name, id)
#file = input_files[1,]
dir.create("out")
dir.create("out/cell_cycle")
dir.create("out/full_dataset")
#write.csv(as.data.frame(processed_files %>% select(name, id)), file="single_cell_runs.csv", quote=T,row.names=F)
for(j in 1:nrow(input_files)){
  #print(file)
  in_folder = input_files$name[j]
  print(in_folder)
  out_folder= glue("out/cell_cycle/{in_folder}")
  out_folder2= glue("out/full_dataset/{in_folder}")
  dir.create(out_folder)
  dir.create(out_folder2)
  matching_results_folder = results_files[results_files$name == in_folder,]
  matching_processed_folder = processed_files[processed_files$name ==in_folder,]
  processed_files_in = drive_ls(matching_processed_folder)
  processed_files_in = processed_files_in[processed_files_in$name == "filtered_feature_bc_matrix",]
  matrix_files = drive_ls(processed_files_in)
  for(i in 1:nrow(matrix_files)){
    #print(i)
    file_name = matrix_files$name[i]
    id = matrix_files$id[i]
    drive_download(matrix_files[i,],path=glue("tmp/{file_name}"),overwrite = T)
  }
  list_folder_processed = drive_ls(matching_results_folder)
  # Check whether diploid of haploid experiment #
  #  m
  cellfilter_rds = "cellFilter.RDS"
  if(sum(list_folder_processed$name == "cell_cycle") == 0){
    cell_cycle_folder = drive_mkdir(path = matching_results_folder, name = "cell_cycle",overwrite = F)
  }else{
    cell_cycle_folder = list_folder_processed[list_folder_processed$name == "cell_cycle",]
  }
  #cell_cycle_contents = drive_ls(cell_cycle_folder)
  #drive_rm(cell_cycle_contents)
  cell_filter = list_folder_processed[list_folder_processed$name == cellfilter_rds,]
  #drive_mkdir(path = matching_results_folder, name = "cell_cycle",overwrite = F)
  drive_download(cell_filter,path=glue("tmp/{cellfilter_rds}"),overwrite = T)
  #drive_get("covid/from_scratch/")
  cell_filter = readRDS(glue("tmp/{cellfilter_rds}"))
  if(file.exists(paste("tmp/","barcordes.tsv.gz",sep=""))){
    cc = cellcycle_annotation("tmp/",cell_filter)
  }else{
    out = (paste("tmp/","barcodes.tsv.gz",sep=""))
    in_f = (paste("tmp/","barcodes.tsv",sep=""))
    system(glue("gzip -c {in_f} > {out}"))
    cc = cellcycle_annotation("tmp/",cell_filter)
  }
  #all = cc$all_list
  #cc = cc$cc 
  #cc_tmp = list()
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
  #drive_upload(out_png, path=cell_cycle_folder,name="clustering.png")
  out_rds_cell_cycle = glue("{out_folder}/cell_cycle.RDS")
  saveRDS(cc,file=out_rds_cell_cycle)
  # out_all_rds_cell_cycle = glue("{out_folder}/all.RDS")
  # saveRDS(all,file=out_all_rds_cell_cycle)
  #drive_upload("tmp/cell_cycle.RDS",path=cell_cycle_folder,name="cell_cycle.RDS",type = "application/zip")
  all_data = read_and_normalize("tmp/",cell_filter)
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

#j = 1 
j = 1
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
DimPlot(cell_rds$cc_seurat,label = T,label.size = 24)
#cell_rds$cc2
#FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G1"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G1"
cell_rds$cc2 %>% filter(cluster == 7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G2/M"
cell_rds$cc_seurat$cell_cycle = cell_cycle
DimPlot(cell_rds$cc_seurat,label = T,label.size = 24,group.by="cell_cycle")
cell_rds$cc_seurat@meta.data %>% ggplot(aes(x=cell_cycle)) + geom_histogram(stat="count")
make_cc_assignment_plots_and_write(cell_rds,name)


j = 2
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
DimPlot(cell_rds$cc_seurat,label = T,label.size = 24)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
#FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G1"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 ==2] = "ALPHA"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G1"
cell_rds$cc2 %>% filter(cluster == 8)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 8] = "G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
DimPlot(cell_rds$cc_seurat,label = T,label.size = 24,group.by="cell_cycle")
cell_rds$cc_seurat@meta.data %>% ggplot(aes(x=cell_cycle)) + geom_histogram(stat="count")
cell_df = data.frame(cell_name=colnames(cell_rds$cc_seurat),cell_cycle=cell_rds$cc_seurat@meta.data$cell_cycle,seurat_clusters=cell_rds$cc_seurat@meta.data$SCT_snn_res.0.3)
make_cc_assignment_plots_and_write(cell_rds,name)

j = 3
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
#FeaturePlot(cell_rds$cc_seurat,anno$Marker)

#DimPlot(cell_rds$cc_seurat,label = T,label.size = 24)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)


j = 4 
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])

all_rds = get_all_rds(input_files$name[j])
FeaturePlot(cell_rds$cc_seurat,anno$Marker)
DimPlot(cell_rds$cc_seurat,label = T,label.size = 24)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G1"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "ALPHA"
cell_rds$cc_seurat$cell_cycle = cell_cycle
DimPlot(cell_rds$cc_seurat,label = T,label.size = 24,group.by="cell_cycle")
make_cc_assignment_plots_and_write(cell_rds,name)


j = 5
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
FeaturePlot(cell_rds$cc_seurat,anno$Marker)
DimPlot(cell_rds$cc_seurat,label = T,label.size = 24)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G1"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "M/G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "M/G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
DimPlot(cell_rds$cc_seurat,label = T,label.size = 24,group.by="cell_cycle")
make_cc_assignment_plots_and_write(cell_rds,name)



j = 6 
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
#FeaturePlot(cell_rds$cc_seurat,anno$Marker)

cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G1"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "S"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G1"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G1/S"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)
### Skip 7 and 8
j = 7
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
DimPlot(cell_rds$cc_seurat)

cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "M/G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
table(cell_cycle)
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)

#print(input_files$name[j])
#name = input_files$name[j]
j = 8
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "M/G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc_seurat$cell_cycle = cell_cycle 
make_cc_assignment_plots_and_write(cell_rds,name)


j = 9
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
#FeaturePlot(cell_rds$cc_seurat,anno$Marker)

cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G1/S"
#cell_rds$cc2 %>% filter(cluster == 6)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G1/S"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)

j = 10
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
#FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G1"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "S"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G1/S"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)
## Weird clusters ##
j = 11
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "S"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 8)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 8] = "G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)

j = 12
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))

cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G2/M"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)


#FeaturePlot(cell_rds$cc_seurat,"FIT3")
j = 13
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
#FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "S"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "ALPHA"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
#cell_rds$cc2 %>% filter(cluster == 4)
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)

j = 14
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
#FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G1"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "ALPHA"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "M/G1"
#cell_rds$cc2 %>% filter(cluster == 7)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "M/G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)


### without 

#cell_rds$assignmennts 

j = 15
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
#FeaturePlot(cell_rds$cc_seurat,"PCL1")

cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
DimPlot(c)
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G1"
cell_rds$cc2 %>% filter(cluster == 7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G1/S"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)
#FeaturePlot(cell_rds$cc_seurat,"HO")

j = 16
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G1"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G1/S"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cell_rds,name)
### Diploids
input_files = read.csv("single_cell_runs.csv")
input_files = input_files %>% filter(type == "diploids" & good == "Y")
folder_in  = drive_find(q = "sharedWithMe and name = 'yeast'",n_max = 25)
results=drive_ls(folder_in)[drive_ls(folder_in)[,1] == "results",]
processed =drive_ls(folder_in)[drive_ls(folder_in)[,1] == "processed",]
results_files = drive_ls(results)
processed_files = drive_ls(processed)
#processed_files %>% select(name, id)
#file = input_files[1,]
dir.create("out")
dir.create("out/cell_cycle")
dir.create("out/full_dataset")

for(j in 1:nrow(input_files)){
  in_folder = input_files$name[j]
  
  print(in_folder)
  out_folder= glue("out/cell_cycle/{in_folder}")
  out_folder2= glue("out/full_dataset/{in_folder}")
  matching_results_folder = results_files[results_files$name == in_folder,]
  matching_processed_folder = processed_files[processed_files$name ==in_folder,]
  processed_files_in = drive_ls(matching_processed_folder)
  processed_files_in = processed_files_in[processed_files_in$name == "filtered_feature_bc_matrix",]
  matrix_files = drive_ls(processed_files_in)
  for(i in 1:nrow(matrix_files)){
    #print(i)
    file_name = matrix_files$name[i]
    id = matrix_files$id[i]
    drive_download(matrix_files[i,],path=glue("tmp/{file_name}"),overwrite = T)
  }
  if(file.exists(paste("tmp/","barcordes.tsv.gz",sep=""))){
    #cc = cellcycle_annotation("tmp/",cell_filter)
    print("Barcodes_exist")
  }else{
    out = (paste("tmp/","barcodes.tsv.gz",sep=""))
    in_f = (paste("tmp/","barcodes.tsv",sep=""))
    system(glue("gzip -c {in_f} > {out}"))
  }
  
  diploid_assignments = "diploid_assignments.RDS"
  list_folder_processed = drive_ls(matching_results_folder)
  cell_filter = list_folder_processed[list_folder_processed$name == diploid_assignments,]
  drive_download(cell_filter,path=glue("tmp/{diploid_assignments}"),overwrite = T)
  diploid_assignments = readRDS(glue("tmp/{diploid_assignments}"))
  diploid_annotation(in_seurat = "tmp/",diploid_assignments = diploid_assignments,in_folder = in_folder)
}

#### samples ####
j = 1
name  = input_files$name[j]
in_dir = glue("out/cell_cycle/{name}/")
diploid_names  = list.files(in_dir,full.names = T)[grep(".png$",list.files(in_dir,full.names = T), invert=T)]
length(diploid_names)

#print(name)
#print(diploid_names)
i = 1

plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
#FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 ==6] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 ==7] = "G2/M"
#cell_rds$cc2 %>% filter(cluster == 8)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 ==8] = "HAPLOIDS"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

i = 2
### Looks like these might be haploid segregants...
cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))



plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

#FeaturePlot(cell_rds$cc_seurat,anno$Marker)
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))

cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "HAPLOIDS"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "HAPLOIDS"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "HAPLOIDS"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "HAPLOIDS"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "HAPLOIDS"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

i = 3 
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "S"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "HAPLOIDS"
cell_rds$cc_seurat$cell_cycle = cell_cycle
#DimPlot(cell_rds$cc_seurat)
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

i = 4
### Weird probably actually a haploid ###
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "HAPLOIDS"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "HAPLOIDS"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "HAPLOIDS"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "HAPLOIDS"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "HAPLOIDS"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "HAPLOIDS"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])
#cell_rds$cc2 %>% filter(cluster == 6)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "HAPLOIDS"

j = 2
name  = input_files$name[j]

print(name)
in_dir = glue("out/cell_cycle/{name}/")
diploid_names  = list.files(in_dir,full.names = T)[grep(".png$",list.files(in_dir,full.names = T), invert=T)]
i = 1
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "S"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "M/G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

i = 2
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "S"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G1/S"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

i = 3 
diploid_names[i]

plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
cell_cycle = rep("HAPLOIDS",ncol(cell_rds$cc_seurat))
### HALPOIDS ####
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

i = 4 
diploid_names[i]
cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))

cell_cycle = rep("HAPLOIDS",ncol(cell_rds$cc_seurat))
### HALPOIDS ####
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])


j  = 3
i = 1


name  = input_files$name[j]
name
in_dir = glue("out/cell_cycle/{name}/")
diploid_names  = list.files(in_dir,full.names = T)[grep(".png$",list.files(in_dir,full.names = T), invert=T)]
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

#i = 1
cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_rds$cc2 %>% filter(cluster ==5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])
#cell_rds$cc2 %>% filter(cluster == 5)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
#cell_rds$cc2 %>% filter(cluster == 6)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 ==6] = "M/G1"
#cell_rds$cc2 %>% filter(cluster == 7)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 ==7] = "G2/M"
#cell_rds$cc2 %>% filter(cluster == 8)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 ==8] = "HAPLOIDS"

i = 2

plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

#name  = input_files$name[j]
#in_dir = glue("out/cell_cycle/{name}/")
#diploid_names  = list.files(in_dir,full.names = T)[grep(".png$",list.files(in_dir,full.names = T), invert=T)]
#i = 1
cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "S"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_rds$cc2 %>% filter(cluster ==5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

i = 3
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "S"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster ==5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster ==6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G2/M"
cell_rds$cc2 %>% filter(cluster ==7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "M/G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

i = 4
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster ==5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster ==6)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "ALPHA"
#cell_rds$cc2 %>% filter(cluster ==7)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G2/M"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

j  = 4
i = 1
name  = input_files$name[j]
in_dir = glue("out/cell_cycle/{name}/")
diploid_names  = list.files(in_dir,full.names = T)[grep(".png$",list.files(in_dir,full.names = T), invert=T)]
cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "S"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G2/M"
cell_rds$cc2 %>% filter(cluster ==5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster ==6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "M/G1"
cell_rds$cc2 %>% filter(cluster ==7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G1/S"
cell_rds$cc2 %>% filter(cluster ==8)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 8] = "G1/S"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

j  = 5
i = 1
#### Looks like there is no G1????
name  = input_files$name[j]
in_dir = glue("out/cell_cycle/{name}/")
diploid_names  = list.files(in_dir,full.names = T)[grep(".png$",list.files(in_dir,full.names = T), invert=T)]
cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])


cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "S"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "S"
cell_rds$cc2 %>% filter(cluster ==5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G1/S"
cell_rds$cc2 %>% filter(cluster ==6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G2/M"
cell_rds$cc2 %>% filter(cluster ==7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G1/S"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

#j  = 5
i = 2
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "S"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster ==5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list =cell_rds, name, folder=diploid_names[i])

i = 3
cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
#cell_rds$cc2 %>% filter(cluster ==5)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list  = cell_rds, name, folder=diploid_names[i])

#cell_rds$cc2 %>% filter(cluster ==6)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G2/M"

j = 6
i = 1

plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

name  = input_files$name[j]
in_dir = glue("out/cell_cycle/{name}/")
diploid_names  = list.files(in_dir,full.names = T)[grep(".png$",list.files(in_dir,full.names = T), invert=T)]
cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "G2/M"
#cell_rds$cc2 %>% filter(cluster == 7)
#cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G2/M"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

#FeaturePlot(cell_rds$cc_seurat,"MFA1")

i = 2
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))

#cell2 = readRDS(glue("{diploid_names[i]}/all.RDS"))
diploid_names[i]

cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "HAPLOIDS"
cell_rds$cc2 %>% filter(cluster == 7)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G2/M"
cell_rds$cc_seurat$cell_cycle = cell_cycle
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle = "HAPLOIDS"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

i = 3
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "S"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "M/G1"
#cell_rds$cc2 %>% filter(cluster == 7)
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])


i = 4
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 2)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "M/G1"
cell_rds$cc2 %>% filter(cluster == 3)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "S"
cell_rds$cc2 %>% filter(cluster == 4)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_rds$cc2 %>% filter(cluster == 5)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 6)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "HAPLOIDS"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])


#cell_rds$cc2 %>% filter(cluster ==5)
j = 7
i = 1
name  = input_files$name[j]
in_dir = glue("out/cell_cycle/{name}/")
diploid_names  = list.files(in_dir,full.names = T)[grep(".png$",list.files(in_dir,full.names = T), invert=T)]

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_rds$cc2 %>% filter(cluster == 0)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_rds$cc2 %>% filter(cluster == 1)
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "G2/M"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 7] = "G2/M"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "G1/S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 6] = "M/G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])

i = 2
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

name  = input_files$name[j]
in_dir = glue("out/cell_cycle/{name}/")
diploid_names  = list.files(in_dir,full.names = T)[grep(".png$",list.files(in_dir,full.names = T), invert=T)]

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]
cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "M/G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "G1/S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "M/G1"
cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])


i = 3
plot_MFA1(diploid_names[i], name, folder=diploid_names[i])

name  = input_files$name[j]
in_dir = glue("out/cell_cycle/{name}/")
diploid_names  = list.files(in_dir,full.names = T)[grep(".png$",list.files(in_dir,full.names = T), invert=T)]

cell_rds =readRDS(glue("{diploid_names[i]}/cell_cycle.RDS"))
diploid_names[i]

cell_cycle = rep(NA,ncol(cell_rds$cc_seurat))
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 0] = "G2/M"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 1] = "S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 2] = "G1/S"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 3] = "M/G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 4] = "M/G1"
cell_cycle[cell_rds$cc_seurat$SCT_snn_res.0.3 == 5] = "G1/S"

cell_rds$cc_seurat$cell_cycle = cell_cycle
make_cc_assignment_plots_and_write(cc_list = cell_rds, name, folder=diploid_names[i])



input_files = input_files %>% filter(type == "diploids" & good == "Y")
folder_in  = drive_find(q = "sharedWithMe and name = 'yeast'",n_max = 25)
results=drive_ls(folder_in)[drive_ls(folder_in)[,1] == "results",]
processed =drive_ls(folder_in)[drive_ls(folder_in)[,1] == "processed",]
results_files = drive_ls(results)
processed_files = drive_ls(processed)
#processed_files %>% select(name, id)
#file = input_files[1,]

