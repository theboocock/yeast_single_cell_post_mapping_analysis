library(tidyverse)
library(glue)
library(googledrive)
library(glmnet)
library(cowplot)
library(Matrix)
library(Seurat)
library(tidyverse)
library(rtracklayer)
set.seed(42)


preprocess_data = function(in_seurat){
  expr_obj = Read10X(in_mm[1])
  expr_obj = CreateSeuratObject(expr_obj,min.cells = 3,min.features = 100)
  expr_obj <- NormalizeData(expr_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  expr_obj <- FindVariableFeatures(expr_obj, selection.method = "vst", nfeatures = 2000, verbose = F)
  expr_obj <- FindNeighbors(expr_obj, dims = 1:12)
  expr_obj <- FindClusters(expr_obj, resolution = 0.5)
  
}

preprocess_data_obj = function(expr_obj,sc_transform=T){
  #expr_obj = Read10X(in_mm[1])
  #expr_obj = CreateSeuratObject(expr_obj,min.cells = 3,min.features = 100)
  #expr_obj <- NormalizeData(expr_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  if(sc_transform){
    expr_obj <- SCTransform(expr_obj)
  }else{
    expr_obj = NormalizeData(expr_obj, normalization.method="LogNormalize",scale.factor=10000, verbose=F)
    expr_obj = ScaleData(expr_obj)
  }
  #expr_obj <- ScaleData(expr_obj)
  expr_obj <- FindVariableFeatures(expr_obj, selection.method = "vst", nfeatures = 2000, verbose = F)
  expr_obj <- RunPCA(expr_obj)
  expr_obj <- FindNeighbors(expr_obj, dims = 1:12)
  expr_obj <- FindClusters(expr_obj, resolution = 0.3)
  expr_obj <- RunUMAP(expr_obj,dims=1:2)
  return(expr_obj)
}  



integrate_data = function(rdss,cell_cycle_annotations,idxs=NULL){
  # Not finished#
  if(!is.null(idxs)){
    rdss = rdss[idxs]
    cell_cycle_annotations =cell_cycle_annotations[idxs]
  }
  for (i in 1:length(rdss)){
    x =readRDS(rdss[i])
    seurat_all_list[[i]] = x$cc_seurat
    a = read.csv(cell_cycle_annotations[i])
    cell_cycle_annotation_list[[i]] = a
  }
  
  
  
  
}

crosses.to.parents=list(
  '375'=c("M22", "BYa"),           #1
  'A'  =c("BYa", "RMx"),           #2
  '376'=c("RMx", "YPS163a"),       #3
  'B'  =c("YPS163a", "YJM145x"),   #4
  '377'=c("YJM145x", "CLIB413a"),  #5
  '393'=c("CLIB413a", "YJM978x"),  #6
  '381'=c("YJM978x", "YJM454a"),   #7
  '3008'=c("YJM454a", "YPS1009x"),  #8
  '2999'=c("YPS1009x", "I14a"),     #9
  '3000'=c("I14a", "Y10x"),         #10
  '3001'=c("Y10x", "PW5a"),         #11
  '3049'=c("PW5a", "273614xa"),     #12
  '3003'=c("273614xa", "YJM981x"),  #13
  '3004'=c("YJM981x", "CBS2888a"),  #14
  '3043'=c("CBS2888a", "CLIB219x"), #15
  '3028'=c("CLIB219x", "M22")       #16
)
chroms=paste0('chr', as.roman(1:16)) 
source("load_cc.R")
read_and_normalize = function(in_seurat,filt_cells,gene_subset=NULL,n_pcs=12){
  expr_obj2= Read10X(in_seurat)
  keep = colnames(expr_obj2) %in% names(filt_cells)[filt_cells]
  if(is.null(gene_subset)){
    expr_obj2 = expr_obj2[,keep]
  }else{
    expr_obj2 = expr_obj2[rownames(expr_obj2) %in% gene_subset,keep]
  }
  expr_obj2 = CreateSeuratObject(expr_obj2)
  expr_obj_sctransform = SCTransform(expr_obj2)
  expr_obj_sctransform <- RunPCA(expr_obj_sctransform)
  #expr_obj_sctransform <- JackStraw(expr_obj_sctransform, num.replicate = 100)
  expr_obj_sctransform <- FindNeighbors(expr_obj_sctransform, dims = 1:n_pcs)
  expr_obj_sctransform <- FindClusters(expr_obj_sctransform, resolution = .3)
  expr_obj_sctransform <- RunUMAP(expr_obj_sctransform, dims = 1:n_pcs)
  markers = FindAllMarkers(expr_obj_sctransform)
  return(list(seurat_obj=expr_obj_sctransform,markers=markers))
}



cellcycle_annotation = function(in_seurat,filt_cells,avg_log2_filt=0.2,p_filt=0.05,n_pcs=12,resolution=0.3){
  print(in_seurat)
  expr_obj2= Read10X(in_seurat)
  keep = colnames(expr_obj2) %in% names(filt_cells)[filt_cells]
  expr_obj2 = expr_obj2[rownames(expr_obj2) %in% cell_cycle_df$NAME,keep]
  expr_obj2 = CreateSeuratObject(expr_obj2)
  expr_obj_sctransform = SCTransform(expr_obj2)
  expr_obj_sctransform <- RunPCA(expr_obj_sctransform)
  #expr_obj_sctransform <- JackStraw(expr_obj_sctransform, num.replicate = 100)
  expr_obj_sctransform <- FindNeighbors(expr_obj_sctransform, dims = 1:n_pcs)
  expr_obj_sctransform <- FindClusters(expr_obj_sctransform, resolution = resolution)
  expr_obj_sctransform <- RunUMAP(expr_obj_sctransform, dims = 1:n_pcs)
  ### Opps make cc2
  markers = FindAllMarkers(expr_obj_sctransform)
  cc2 = markers %>% filter((avg_log2FC) > avg_log2_filt &p_val_adj < p_filt) %>% inner_join(cell_cycle_big_df,by=c("gene"="NAME")) 
  annotate = cc2 %>% group_by(cluster) %>% summarise(cc=names(sort(table(Peak), decreasing = TRUE))[1])
  cell_cycle = rep(NA,length(expr_obj_sctransform$SCT_snn_res.0.3))
  i = 1
  for(clust in annotate$cluster){
    cell_cycle[expr_obj_sctransform$SCT_snn_res.0.3 == clust] = annotate$cc[i]
    i = i +1
  }
  expr_obj_sctransform@meta.data$cell_cycle = cell_cycle
  cc_df = data.frame(cell_id=rownames(expr_obj_sctransform@meta.data),
                     cell_cycle=expr_obj_sctransform@meta.data$cell_cycle,
                     unsupervised=expr_obj_sctransform@meta.data$SCT_snn_res.0.3)
  return(list(cc_seurat=expr_obj_sctransform,assignmennts=cell_cycle, cc_df=cc_df,cc2=cc2))
  # all_list = list(all_seurat=expr_obj_sctransform_all,all_markers=all_markers)
  #return(list(all_list=all_list, cc=cc)) 
}
avg_log2_filt=0.2
p_filt=0.05
n_pcs=12
resolution=0.3
likelihood_diff_filter=200

diploid_annotation = function(in_seurat,diploid_assignments,in_folder,haploid=F,avg_log2_filt=0.2,p_filt=0.05,n_pcs=12,resolution=0.3,likelihood_diff_filter=200){
  #print(in_seurat)
  expr_obj2= Read10X(in_seurat)
  #print("HEH")
  keep = diploid_assignments %>% filter(is.na(diploid_assignment_likdiff) | diploid_assignment_likdiff > likelihood_diff_filter) %>% select(barcode)
  diploid_assignments_filt  = diploid_assignments %>% filter(is.na(diploid_assignment_likdiff)  | diploid_assignment_likdiff > likelihood_diff_filter)
  #keep = colnames(expr_obj2) %in% names(filt_cells)[filt_cells]
  print(head(diploid_assignments_filt))
  expr_obj3 = expr_obj2[,keep$barcode]
  
  expr_obj2 = expr_obj2[rownames(expr_obj2) %in% cell_cycle_df$NAME,keep$barcode]
  if(!haploid){
    expr_obj2 = expr_obj2[!rownames(expr_obj2) %in% str_replace_all(scer_hsg_genes_cell_cycle,"SCER-",""),]
  }
  
  # # # # # # # #
  seurat_objects = list()
  plots = list()
  out_dir = glue("out/cell_cycle/{in_folder}")
  dir.create(out_dir)
  i = 1
  print("two")
  for(diploid_tmp_name in unique(diploid_assignments_filt$diploid_name)){
    print(diploid_tmp_name)
    out_dir_tmp = glue("{out_dir}/{diploid_tmp_name}")
    dir.create(out_dir_tmp)
    #out_dir_tmp = glue("out/cell_cycle/{in_folder}")
    keep_diploid = diploid_assignments_filt %>% dplyr::filter(diploid_name == diploid_tmp_name)
    expr_obj_tmp = expr_obj2[rownames(expr_obj2) %in% cell_cycle_df$NAME,keep_diploid$barcode]
    expr_obj_tmp = CreateSeuratObject(expr_obj_tmp)
    expr_obj_sctransform = SCTransform(expr_obj_tmp)
    expr_obj_sctransform <- RunPCA(expr_obj_sctransform)
    #expr_obj_sctransform <- JackStraw(expr_obj_sctransform, num.replicate = 100)
    expr_obj_sctransform <- FindNeighbors(expr_obj_sctransform, dims = 1:n_pcs)
    expr_obj_sctransform <- FindClusters(expr_obj_sctransform, resolution = resolution)
    expr_obj_sctransform <- RunUMAP(expr_obj_sctransform, dims = 1:n_pcs)
    DimPlot(expr_obj_sctransform)
    
    
    expr_obj_tmp2 = expr_obj3[,keep_diploid$barcode]
    expr_obj_tmp2 = CreateSeuratObject(expr_obj_tmp2)
    expr_obj_sctransform2 = SCTransform(expr_obj_tmp2)
    expr_obj_sctransform2<- RunPCA(expr_obj_sctransform2)
    #expr_obj_sctransform <- JackStraw(expr_obj_sctransform, num.replicate = 100)
    expr_obj_sctransform2 <- FindNeighbors(expr_obj_sctransform2, dims = 1:n_pcs)
    expr_obj_sctransform2 <- FindClusters(expr_obj_sctransform2, resolution = resolution)
    expr_obj_sctransform2 <- RunUMAP(expr_obj_sctransform2, dims = 1:n_pcs)
    markers = FindAllMarkers(expr_obj_sctransform)
    all_out = list(markers=markers, seurat_obj=expr_obj_sctransform2)
    out_rds = glue("{out_dir_tmp}/all.RDS")
    saveRDS(all_out, file=out_rds)
    
    png_dim = glue("{out_dir_tmp}/dimplot.png")
    png(png_dim, width=800, height=600)
    print(DimPlot(expr_obj_sctransform,label=T,label.size = 24))
    dev.off()
    
    plots[[i]] = (DimPlot(expr_obj_sctransform))
    png_cluster = glue("{out_dir_tmp}/feature_plot_cc.png")
    png(png_cluster,width=1600, height=1200)
    print(FeaturePlot(expr_obj_sctransform,anno$Marker))
    dev.off()
    markers = FindAllMarkers(expr_obj_sctransform)
    if(nrow(markers) != 0){
      cc2 = markers %>% filter((avg_log2FC) > avg_log2_filt &p_val_adj < p_filt) %>% inner_join(cell_cycle_big_df,by=c("gene"="NAME")) 
    }else{
      cc2 = markers
    }
    #annotate = cc2 %>% group_by(cluster) %>% summarise(cc=names(sort(table(Peak), decreasing = TRUE))[1])
    #cell_cycle = rep(NA,length(expr_obj_sctransform$SCT_snn_res.0.3))
    #markers = FindAllMarkers(expr_obj_sctransform)
    #cc2 = markers %>% filter((avg_log2FC) > avg_log2_filt &p_val_adj < p_filt) %>% inner_join(cell_cycle_big_df,by=c("gene"="NAME")) 
    cc_out = list(cc_seurat=expr_obj_sctransform,cc2=cc2)
    
    #full_object = 
    #annotate = cc2 %>% group_by(cluster) %>% summarise(cc=names(sort(table(Peak), decreasing = TRUE))[1])
    #cell_cycle = rep(NA,length(expr_obj_sctransform$SCT_snn_res.0.3))
    out_rds = glue("{out_dir_tmp}/cell_cycle.RDS")
    saveRDS(cc_out, file=out_rds)
    
    ## ###
    
    
    
    i = i + 1
  }
  png_one = glue("{out_dir}/all_cluster.png")
  png(png_one,width=1600, height=1200)
  print(cowplot::plot_grid(plotlist = plots,labels = unique(diploid_assignments_filt$diploid_name)))
  dev.off()
  #return(list(cc_seurat=expr_obj_sctransform,assignmennts=cell_cycle, cc_df=cc_df,cc2=cc2)) 
}


#run_multinom_regression

#cellcycle_reannotate = function(in_seurat, avg_log2_filt=0.2,p_filt=0.05){
#  in_seurat@cc2 %>% filter(avg_log2FC > avg_log2_filt) %>% 
#}
get_cell_cycle_rds= function(name){
  return(readRDS(glue("out/cell_cycle/{name}/cell_cycle.RDS")))
}
get_all_rds = function(name){
  return(readRDS(glue("out/cell_cycle/{name}/all_seurat.RDS")))
}
get_cell_cycle_annotations= function(name){
  return(read.csv(glue("out/cell_cycle/{name}/cell_cycle_assignments.csv")))
}
add_cell_cycle_annotation_to_seurat_objects_and_load = function(name){
  cell_cycle_seurat = readRDS(glue("out/cell_cycle/{name}/cell_cycle.RDS"))
  all_rds = readRDS(glue("out/cell_cycle/{name}/all_seurat.RDS"))
  cell_cycle = read.csv(glue("out/cell_cycle/{name}/cell_cycle_assignments.csv"))
  cell_cycle_seurat$cc_seurat$cell_cycle_final = cell_cycle$cell_cycle  
  all_rds$seurat_obj$cell_cycle_final = cell_cycle$cell_cycle
  return(list(all_seurat=all_rds$seurat_obj,cell_cycle_seurat=cell_cycle_seurat$cc_seurat))
}


make_cc_assignment_plots_and_write_ms= function(cc_list, out_folder ){
  cc_object = cc_list$cc_seurat
  png_one=glue("{out_folder}/side_by_side.png") 
  print(png_one)
  png(png_one,width = 1600, height=1200)
  p1 = DimPlot(cc_object,label = T,label.size = 24)
  p2 = DimPlot(cc_object,label = T,label.size = 12, group.by = "cell_cycle")
  print(cowplot::plot_grid(p1,p2))
  dev.off()
  png_two=glue("{out_folder}/cell_cycle_histogram.png") 
  png(png_two, width = 800,height=600)
  print(cc_object@meta.data %>% ggplot(aes(x=cell_cycle)) + geom_histogram(stat="count"))
  dev.off()  
  cell_df = data.frame(cell_name=colnames(cc_object),cell_cycle=cc_object@meta.data$cell_cycle,seurat_clusters=cc_object@meta.data$SCT_snn_res.0.3)
  out_csv = glue("{out_folder}/cell_cycle_assignments.csv")
  write.csv(cell_df, file=out_csv,quote=T)    
  png_rna_count = glue("{out_folder}/ncount_rna.png")
  png(png_rna_count, width=800, height=600)
  print(FeaturePlot(cc_object,features = "nCount_RNA"))
  dev.off()
  out_cc_rds = glue("{out_folder}/cell_cycle_final.RDS")
  saveRDS(cc_object)
  
}
make_cc_assignment_plots_and_write= function(cc_list, name,folder=NULL){
  cc_object = cc_list$cc_seurat
  if (!is.null(folder)){
    out_dir = glue("{folder}/")
  }else{
    out_dir=glue("out/cell_cycle/{name}/")
  }
  png_one=glue("{out_dir}side_by_side.png") 
  print(png_one)
  png(png_one,width = 1600, height=1200)
  p1 = DimPlot(cc_object,label = T,label.size = 24)
  p2 = DimPlot(cc_object,label = T,label.size = 12, group.by = "cell_cycle")
  print(cowplot::plot_grid(p1,p2))
  dev.off()
  png_two=glue("{out_dir}cell_cycle_histogram.png") 
  png(png_two, width = 800,height=600)
  print(cc_object@meta.data %>% ggplot(aes(x=cell_cycle)) + geom_histogram(stat="count"))
  dev.off()  
  cell_df = data.frame(cell_name=colnames(cc_object),cell_cycle=cc_object@meta.data$cell_cycle,seurat_clusters=cc_object@meta.data$SCT_snn_res.0.3)
  out_csv = glue("{out_dir}cell_cycle_assignments.csv")
  write.csv(cell_df, file=out_csv,quote=T)    
  png_three = glue("{out_dir}cell_cycle_feature_plot.png")
  print(png_three)
  png(png_three, width=1600, height=1000)
  print(FeaturePlot(cc_object,anno$Marker))
  dev.off()
  png_rna_count = glue("{out_dir}ncount_rna.png")
  png(png_rna_count, width=800, height=600)
  print(FeaturePlot(cc_object,features = "nCount_RNA"))
  dev.off()  
  out_cc_rds = glue("{out_dir}cell_cycle_final.RDS")
  print(out_cc_rds)
  saveRDS(cc_list, file=out_cc_rds)
  
  
}
plot_MFA1 = function(diploid_names,name, folder=NULL){
  if (!is.null(folder)){
    out_dir = glue("{folder}/")
  }else{
    out_dir=glue("out/cell_cycle/{name}/")
  }
  cell2 = readRDS(glue("{diploid_names}/all.RDS")) 
  png(glue("{out_dir}/mating_factor.png"),width=1600,height=1200)
  print(FeaturePlot(cell2$seurat_obj,str_replace_all(scer_hsg_genes_cell_cycle,"SCER-","")))
  dev.off()
}


get_data_dirs = function(x){
  folder_in  = drive_find(q = "sharedWithMe and name = 'yeast'",n_max = 25)
  results=drive_ls(folder_in)[drive_ls(folder_in)[,1] == "results",]
  processed =drive_ls(folder_in)[drive_ls(folder_in)[,1] == "processed",]
  results_files = drive_ls(results)
  processed_files = drive_ls(processed)
  return(list(results=results_files, processed=processed_files))
}

get_processed_dir = function(x){
  folder_in  = drive_find(q = "sharedWithMe and name = 'yeast'",n_max = 25)
  results=drive_ls(folder_in)[drive_ls(folder_in)[,1] == "results",]
  processed =drive_ls(folder_in)[drive_ls(folder_in)[,1] == "processed",]
  results_files = drive_ls(results)
  return(results_files)
}

