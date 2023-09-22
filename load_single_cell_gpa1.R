library(Seurat)
set.seed(42)
sample_a = Read10X("data/single_cell/sample_a/outs//filtered_feature_bc_matrix//")
sample_b = Read10X("data/single_cell/sample_b//outs/filtered_feature_bc_matrix/")
keep_genes= intersect(rownames(sample_a),rownames(sample_b))
idx_keep_a  = (colSums(sample_a) < 25e3)
idx_keep_b = (colSums(sample_b) < 25e3)
sample_a = sample_a[,idx_keep_a]
sample_b= sample_b[,idx_keep_b]
sample_a_m = CreateSeuratObject(counts=sample_a)
sample_b_m = CreateSeuratObject(counts=sample_b)
gpa_m = merge(sample_a_m,sample_b_m,add.cell.ids=c("A","B"),project="GPA1")
gpa_m = SCTransform(gpa_m)
gpa_m <- RunPCA(gpa_m, verbose = FALSE)
gpa_m <- RunUMAP(gpa_m, dims = 1:30, verbose = FALSE)
gpa_m <- FindNeighbors(gpa_m, dims = 1:30, verbose = FALSE)
gpa_m <- FindClusters(gpa_m, verbose = FALSE)
library(stringr)
sample = unlist(lapply(str_split(rownames(gpa_m@meta.data),"_"), function(x){x[1]}))
gpa_m$sample = sample
#aaa = FindMarkers(gpa_m,group.by="sample",ident.1="A",ident.2 = "B")
#aaa[rownames(aaa) == "MFA1",]
mk = FindMarkers(gpa_m,group.by="sample",ident.1="A",ident.2 = "B",logfc.threshold = 0,min.pct = 0)
#pm = cross_data$B$trans$hotspot_enrichments_and_overlaps$`chrVIII:46887-140660_51`$peak_merged
#pm$gene_name[!pm$gene_name %in% rownames(mk)]
#bl = mk[rownames(mk) %in% pm$gene_name,]  %>% rownames_to_column(var="gene") %>% inner_join(pm,by=c("gene"="gene_name"))
library(ggpubr)
#bl %>% ggplot(aes(y=avg_log2FC,x=Beta.y)) + geom_point() + stat_cor(method="spearman")
sample_a_cc = sample_a[rownames(sample_a) %in% cell_cycle_df$NAME,]
sample_b_cc = sample_b[rownames(sample_b) %in% cell_cycle_df$NAME,]
sample_a_c = CreateSeuratObject(counts=sample_a_cc)
sample_b_c = CreateSeuratObject(counts=sample_b_cc)
gpa_c = merge(sample_a_c,sample_b_c,add.cell.ids=c("A","B"),project="GPA1")
gpa_c = SCTransform(gpa_c)
gpa_c <- RunPCA(gpa_c, verbose = FALSE)
gpa_c <- RunUMAP(gpa_c, dims = 1:30, verbose = FALSE)
gpa_c <- FindNeighbors(gpa_c, dims = 1:12, verbose = FALSE)
gpa_c <- FindClusters(gpa_c, verbose = FALSE,resolution = 0.3)
gpa_c$dataset = unlist(lapply(str_split(rownames(gpa_c@meta.data),"_"), function(x){x[1]}))
gpa_mk = FindAllMarkers(gpa_c)
gpa_mk_df = gpa_mk %>% as_data_frame()
cell_cycle_clusters = rep(NA,length(gpa_c$seurat_clusters))

gpa_mk_df %>% filter(cluster == 0)
cell_cycle_clusters[gpa_c$seurat_clusters ==0]= "S"
gpa_mk_df %>% filter(cluster == 1)
cell_cycle_clusters[gpa_c$seurat_clusters ==1] = "G1"
gpa_mk_df %>% filter(cluster == 2)
cell_cycle_clusters[gpa_c$seurat_clusters ==2] = "G2/M"
gpa_mk_df %>% filter(cluster == 3)
cell_cycle_clusters[gpa_c$seurat_clusters ==3]="G2/M"
gpa_mk_df %>% filter(cluster == 4)
cell_cycle_clusters[gpa_c$seurat_clusters ==4]="M/G1"
gpa_mk_df %>% filter(cluster == 5)
cell_cycle_clusters[gpa_c$seurat_clusters ==5] = "M/G1"
gpa_mk_df %>% filter(cluster == 6)
cell_cycle_clusters[gpa_c$seurat_clusters ==6] = "S"
gpa_mk_df %>% filter(cluster == 7)
cell_cycle_clusters[gpa_c$seurat_clusters ==7] = "G1"
gpa_mk_df %>% filter(cluster == 8)
cell_cycle_clusters[gpa_c$seurat_clusters ==8] = "G1/S"
gpa_mk_df %>% filter(cluster == 9)
cell_cycle_clusters[gpa_c$seurat_clusters ==9] = "G1/S"
gpa_mk_df %>% filter(cluster == 10)
cell_cycle_clusters[gpa_c$seurat_clusters ==10] = "G1/S"
gpa_mk_df %>% filter(cluster == 11)
cell_cycle_clusters[gpa_c$seurat_clusters ==11] = "G2/M"
cell_cycle_clusters[gpa_c$seurat_clusters ==12] = "S"
#table(cell_cycle_clusters)
gpa_c$cell_cycle = cell_cycle_clusters

cc_mfa = cell_cycle_df[cell_cycle_df$NAME != "MFA1",]
sample_a_cc = sample_a[rownames(sample_a) %in% cc_mfa$NAME,]
sample_b_cc = sample_b[rownames(sample_b) %in% cc_mfa$NAME,]
sample_a_cm = CreateSeuratObject(counts=sample_a_cc)
sample_b_cm = CreateSeuratObject(counts=sample_b_cc)
gpa_cm = merge(sample_a_cm,sample_b_cm,add.cell.ids=c("A","B"),project="GPA1")
gpa_cm = SCTransform(gpa_cm)
gpa_cm <- RunPCA(gpa_cm, verbose = FALSE)
gpa_cm <- RunUMAP(gpa_cm, dims = 1:30, verbose = FALSE)
gpa_cm<- FindNeighbors(gpa_cm, dims = 1:12, verbose = FALSE)
gpa_cm <- FindClusters(gpa_cm, verbose = FALSE,resolution = 0.3)
#FeaturePlot(gpa_cm,"MFA1")
gpa_cm$dataset = unlist(lapply(str_split(rownames(gpa_cm@meta.data),"_"), function(x){x[1]}))
gpa_cm@assays = gpa_c@assays
mk_no_mfa1 = FindAllMarkers(gpa_cm)

# Create list of all the important stuff 
cell_cycle_rds_list = list(gpa_cm=gpa_cm, mk_no_mfa1, gpa_c=gpa_c, mk_c_w_mfa1=gpa_mk, gpa_m=gpa_m, mk=mk)




saveRDS(cell_cycle_rds_list,file="../rproj/out/cell_cycle/cell_cycle_plot_object.RDS")



