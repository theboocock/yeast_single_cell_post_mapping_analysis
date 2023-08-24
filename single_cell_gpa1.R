library(Seurat)
set.seed(42)
sample_a = Read10X("data/single_cell/sample_a/outs//filtered_feature_bc_matrix//")
sample_b = Read10X("data/single_cell/sample_b//outs/filtered_feature_bc_matrix/")

keep_genes= intersect(rownames(sample_a),rownames(sample_b))

#sample_a=sample_a[keep_genes,]
#sample_b=sample_b[keep_genes,]

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
#DimPlot(gpa_m,group.by="sample")
aaa = FindMarkers(gpa_m,group.by="sample",ident.1="A",ident.2 = "B")
#aaa[rownames(aaa) == "MFA1",]
mk = FindMarkers(gpa_m,group.by="sample",ident.1="A",ident.2 = "B",logfc.threshold = 0,min.pct = 0)


pm = cross_data$B$trans$hotspot_enrichments_and_overlaps$`chrVIII:46887-140660_51`$peak_merged
#genes_name_trans 

pm$gene_name[!pm$gene_name %in% rownames(mk)]

bl = mk[rownames(mk) %in% pm$gene_name,]  %>% rownames_to_column(var="gene") %>% inner_join(pm,by=c("gene"="gene_name"))
#%>% inner_join(pm,by=c(""))
library(ggpubr)
bl %>% ggplot(aes(y=avg_log2FC,x=Beta.y)) + geom_point() + stat_cor(method="spearman")

View(bl)
mk %>% filt
# A  =420/421
# B = 416/417
fisher.test(base::table(bl$avg_log2FC< 0,bl$Beta.y < 0))

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

p1 = DimPlot(gpa_c, group.by="dataset")
p2=FeaturePlot(gpa_c,c("MFA1"))

cowplot::plot_grid(p1,p2)

gpa_mk = FindAllMarkers(gpa_c)
hist(gpa_c$nCount_RNA)
hist(gpa_m$nCount_RNA,breaks=1000)

### Cell_cycle assignement ###
p1 = DimPlot(gpa_c,label=T,label.size = 12)
p2 = FeaturePlot(gpa_c,anno$Marker)

cowplot::plot_grid(p1,p2)
ggsave("cell_cycle_gpa1.png")

gpa_mk_df = gpa_mk %>% inner_join(cell_cycle_big_df,by=c("gene"="NAME")) %>% filter(avg_log2FC > 1)# %>%
#  group_by(cluster) %>% summarise(table(Peak))


cell_cycle_clusters = rep(NA,length(gpa_c$seurat_clusters))# %>% 

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
table(cell_cycle_clusters)

gpa_c$cell_cycle = cell_cycle_clusters
DimPlot(gpa_c,group.by="cell_cycle")


aaa = table(gpa_c$cell_cycle == "G1",gpa_c$dataset)
aaa
model = list()
out_df_gpa = data.frame()
out_df_gpa1_without_cc = data.frame()
for(cc in unique(gpa_c$cell_cycle)){
  print(cc)
  
  y = as.numeric(gpa_c$cell_cycle == cc)#~ gpa_c$dataset)
  
  aa = glm(y ~ log(gpa_c$nCount_RNA) + gpa_c$dataset,family="binomial")
  
  af = (drop1(aa,test="Chisq"))
  aa1 = glm(y ~ gpa_c$dataset,family="binomial")
  
  
  
  print(summary(aa))
  aa = emmeans::emmeans(aa,spec=~dataset)
  p1 = pairs(aa)
  
  
  out_df_gpa = rbind(out_df_gpa,data.frame(p1) %>% mutate(lrt=af$`Pr(>Chi)`[3],cell_cycle=cc))
  af = (drop1(aa1,test="Chisq"))
  
  aa1 = emmeans::emmeans(aa1, spec=~dataset)
  
  
  p1=pairs(aa1)
  out_df_gpa1_without_cc = rbind(out_df_gpa1_without_cc,data.frame(p1) %>% mutate(lrt=af$`Pr(>Chi)`[2],cell_cycle=cc))
}

out_df_gpa
out_df_gpa1_without_cc




glc = gpa_c$cell_cycle == "G1" | gpa_c$cell_cycle == "M/G1"
aaaa = emmeans::emmeans(glm(glc ~ log(gpa_c$nCount_RNA) +  gpa_c$dataset,family="binomial"),spec=~dataset)

m1  = glm(glc ~ (log(gpa_c$nCount_RNA))+  gpa_c$dataset,family="binomial")


summary(aaaa)

FeaturePlot(gpa_c,"nCount_RNA")


pairs(aaaa)
DimPlot(gpa_c,group.by="dataset")


gpa_c@meta.data %>% group_by(cell_cycle, dataset) %>% summarise(n=n())
#gpa_c$

fisher.test

out_df_gpa %>% ggplot(aes(y=-estimate,x=cell_cycle)) + geom_bar(stat="identity")

apply(aaa,2,sum)

aaa


aaa = table(gpa_c$cell_cycle == "G1",gpa_c$dataset)
fisher.test(aaa)
aaa[2,1]/sum(aaa[,1])
aaa[2,2]/sum(aaa[,2])

1995/(9555 + 1995)
2297/(2297 + 12589)

fisher.test(aaa)


0.1727273/(0.1543061 + 0.1727273)


9555/(9555 + 12589)
1995/(1995 + 2297)



mating_genes = c("AGA1","AGA2","MFA1","STE2","FUS3","PRM5")


bl[bl$gene %in% mating_genes,]

bl %>% filter(gene %in% mating_genes) %>% ggplot(aes(y=-avg_log2FC,x=Beta.y)) + geom_point()


cc_mfa = cell_cycle_df[cell_cycle_df$NAME != "MFA1",]
sample_a_cc = sample_a[rownames(sample_a) %in% cc_mfa$NAME,]
sample_b_cc = sample_b[rownames(sample_b) %in% cc_mfa$NAME,]


2794/(2794 +3248)
8901/(8901 +11916)


sample_a_cm = CreateSeuratObject(counts=sample_a_cc)
sample_b_cm = CreateSeuratObject(counts=sample_b_cc)
gpa_cm = merge(sample_a_cm,sample_b_cm,add.cell.ids=c("A","B"),project="GPA1")


gpa_cm = SCTransform(gpa_cm)
gpa_cm <- RunPCA(gpa_cm, verbose = FALSE)
gpa_cm <- RunUMAP(gpa_cm, dims = 1:30, verbose = FALSE)
gpa_cm<- FindNeighbors(gpa_cm, dims = 1:12, verbose = FALSE)
gpa_cm <- FindClusters(gpa_cm, verbose = FALSE,resolution = 0.3)
FeaturePlot(gpa_cm,"MFA1")
gpa_cm$dataset = unlist(lapply(str_split(rownames(gpa_cm@meta.data),"_"), function(x){x[1]}))


DimPlot(gpa_cm,group.by="dataset",label=T)

DimPlot(gpa_cm,label=T)


aa = table(gpa_cm$dataset,gpa_cm$seurat_clusters)
apply(aa,1,sum)
aa[1,]/11695
aa[2,]/15164


gpa_cm@assays = gpa_c@assays

#all_m = FindAllMarkers(gpa_cm)
all_m %>% filter(cluster == 10) %>% filter(avg_log2FC >1)
out_df_gpa2= data.frame()
for(cc in unique(gpa_cm$seurat_clusters)){
  print(cc)
  aa = glm(gpa_cm$seurat_clusters == cc ~ gpa_cm$dataset,family="binomial")
  print(summary(aa))
  aa = emmeans::emmeans(aa,spec=~dataset)
  p1 = pairs(aa)
  out_df_gpa2 = rbind(out_df_gpa2,data.frame(p1) %>% mutate(cell_cycle=cc))# %>% mutate(estimate=-estimate))
}
summary(glm(gpa_cm$seurat_clusters == 0 ~ gpa_cm$dataset,family = "binomial"))


FeaturePlot(gpa_cm,anno$Marker)
