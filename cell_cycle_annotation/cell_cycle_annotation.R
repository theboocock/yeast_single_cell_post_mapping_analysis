source("pipeline/load_cc.R")
#source("pipeline/utils.R")
source("cell_cycle_annotation/cell_cycle_annotation_fx.R")
# Let's start with the 480 BYxRM #
library(tidyverse)
library(glue)
library(googledrive)
library(cowplot)
# Note: Josh filters cells with > 10,000 UMIs. 
#in_mm = Sys.glob("/HOFFMAN/old_alignments/processed/processed/*/*outs/filtered_feature_bc_matrix/")
input_files = read.csv("data/single_cell_runs.csv")
input_files = input_files %>% filter(type == "haploids")

#j = 1 
j = 1
print(input_files$name[j])
name = input_files$name[j]
cell_rds = get_cell_cycle_rds(input_files$name[j])
DimPlot(cell_rds$cc_seurat,label = T,label.size = 24)
#cell_rds$cc2
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
input_files = read.csv("data/single_cell_runs.csv")
input_files = input_files %>% filter(type == "diploids" & good == "Y")
#file = input_files[1,]
#dir.create("out")
#dir.create("out/cell_cycle")
#dir.create("out/full_dataset")
dir.create("data/cell_cycle/out",recursive = T)

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

