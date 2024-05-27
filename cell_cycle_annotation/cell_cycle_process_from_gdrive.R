source("pipeline/load_cc.R")
#source("pipeline/utils.R")
source("cell_cycle_annotation/cell_cycle_annotation_fx.R")
# Let's start with the 480 BYxRM #
library(tidyverse)
library(glue)
library(googledrive)
library(cowplot)
# Haploids
input_files = read.csv("data/single_cell_runs.csv")
input_files = input_files %>% filter(type == "haploids")
folder_in  = drive_find(q = "sharedWithMe and name = 'yeast'",n_max = 25)
results=drive_ls(folder_in)[drive_ls(folder_in)[,1] == "results",]
processed =drive_ls(folder_in)[drive_ls(folder_in)[,1] == "processed",]
results_files = drive_ls(results)
processed_files = drive_ls(processed)
#processed_files %>% select(name, id)
#file = input_files[1,]
dir.create("data/cell_cycle/out",recursive = T)
dir.create("data/cell_cycle/out/cell_cycle",recursive = T)
dir.create("data/cell_cycle/out/full_dataset",recursive = T)
set.seed(42)
#write.csv(as.data.frame(processed_files %>% select(name, id)), file="single_cell_runs.csv", quote=T,row.names=F)
for(j in 1:nrow(input_files)){
  #print(file)
  in_folder = input_files$name[j]
  print(in_folder)
  out_folder= glue("data/cell_cycle/out/{in_folder}")
  out_folder2= glue("data/cell_cycle/out/full_dataset/{in_folder}")
  dir.create(out_folder)
  #dir.create(out_folder2)
  matching_results_folder = results_files[results_files$name == in_folder,]
  matching_processed_folder = processed_files[processed_files$name ==in_folder,]
  processed_files_in = drive_ls(matching_processed_folder)
  processed_files_in = processed_files_in[processed_files_in$name == "filtered_feature_bc_matrix",]
  matrix_files = drive_ls(processed_files_in)
  for(i in 1:nrow(matrix_files)){
    #print(i)
    file_name = matrix_files$name[i]
    id = matrix_files$id[i]
    drive_download(matrix_files[i,],path=glue("{out_folder}/{file_name}"),overwrite = T)
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
  drive_download(cell_filter,path=glue("{out_folder}/{cellfilter_rds}"),overwrite = T)
  #drive_get("covid/from_scratch/")
  cell_filter = readRDS(glue("{out_folder}/{cellfilter_rds}"))
  if(file.exists(paste("{out_folder}/","barcordes.tsv.gz",sep=""))){
    cc = cellcycle_annotation(glue("{out_folder}/"),cell_filter)
  }else{
    out = (paste(glue("{out_folder}/"),"barcodes.tsv.gz",sep=""))
    in_f = (paste(glue("{out_folder}/"),"barcodes.tsv",sep=""))
    system(glue("gzip -c {in_f} > {out}"))
    cc = cellcycle_annotation(glue("{out_folder}/"),cell_filter)
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
folder_in  = drive_find(q = "sharedWithMe and name = 'yeast'",n_max = 25)
results=drive_ls(folder_in)[drive_ls(folder_in)[,1] == "results",]
processed =drive_ls(folder_in)[drive_ls(folder_in)[,1] == "processed",]
results_files = drive_ls(results)
processed_files = drive_ls(processed)
dir.create("data/cell_cycle/out",recursive = T)
dir.create("data/cell_cycle/out/cell_cycle",recursive = T)
dir.create("data/cell_cycle/out/full_dataset",recursive = T)
for(j in 1:nrow(input_files)){
  in_folder = input_files$name[j]
  out_folder= glue("data/cell_cycle/out/{in_folder}")
  out_folder2= glue("data/cell_cycle/out/full_dataset/{in_folder}")
  matching_results_folder = results_files[results_files$name == in_folder,]
  matching_processed_folder = processed_files[processed_files$name ==in_folder,]
  processed_files_in = drive_ls(matching_processed_folder)
  processed_files_in = processed_files_in[processed_files_in$name == "filtered_feature_bc_matrix",]
  matrix_files = drive_ls(processed_files_in)
  dir.create(out_folder)
  #dir.create(out_folder2)
  for(i in 1:nrow(matrix_files)){
    #print(i)
    file_name = matrix_files$name[i]
    id = matrix_files$id[i]
    drive_download(matrix_files[i,],path=glue("{out_folder}/{file_name}"),overwrite = T)
  }
  if(file.exists(paste(glue("{out_folder}/"),"barcordes.tsv.gz",sep=""))){
    #cc = cellcycle_annotation("tmp/",cell_filter)
    print("Barcodes_exist")
  }else{
    out = (paste(glue("{out_folder}/"),"barcodes.tsv.gz",sep=""))
    in_f = (paste(glue("{out_folder}/"),"barcodes.tsv",sep=""))
    system(glue("gzip -c {in_f} > {out}"))
  }
  
  diploid_assignments = "diploid_assignments.RDS"
  list_folder_processed = drive_ls(matching_results_folder)
  cell_filter = list_folder_processed[list_folder_processed$name == diploid_assignments,]
  drive_download(cell_filter,path=glue("{out_folder}/{diploid_assignments}"),overwrite = T)
  diploid_assignments = readRDS(glue("{out_folder}/{diploid_assignments}"))
  in_seurat = glue("{out_folder}/")
  diploid_assignments = diploid_assignments
  in_folder = in_folder 
  diploid_annotation(in_seurat = glue("{out_folder}/"),diploid_assignments = diploid_assignments,in_folder = in_folder)
}