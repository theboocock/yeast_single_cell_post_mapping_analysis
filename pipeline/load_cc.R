
set.seed(42)
anno = readRDS("data/cc/ccdf.20210511.RDS")
cell_cycle_df = read.csv("data/cc/botstein_genes.csv")

gff_in = import.gff("data/saccharomyces_cerevisiae.gff")

#genes_gtf = import.gff("/media/hoffman//sceqtl/2021/ref/yeast/genes/genes.gtf")
#genes_gtf_df = genes_gtf %>% as_data_frame()
genes_gtf_df = gff_in %>% filter(type=="gene") %>% as_data_frame()
genes_gtf_df$gene[is.na(genes_gtf_df$gene)] = genes_gtf_df$Name[is.na(genes_gtf_df$gene)] 


cell_cycle_df$NAME = genes_gtf_df$gene[match(cell_cycle_df$ORF,genes_gtf_df$Name)]
cell_cycle_df = cell_cycle_df %>% filter(!is.na(NAME))
#cell_cycle_df$NAME = genes_gtf_df$Name[match(cell_cycle_df$ORF,genes_gtf_df$Name)]
genes_name_trans = genes_gtf_df %>% dplyr::select(Name,gene)
colnames(genes_name_trans) = c("gene_id","gene_name")
genes_name_trans$gene_name[is.na(genes_name_trans$gene_name)] = genes_name_trans$gene_id[is.na(genes_name_trans$gene_name)]
new_names = genes_gtf_df$gene[match(anno$Marker, genes_gtf_df$Name)]
anno$Marker = new_names
anno$Type[anno$Type == "G1"] = "G1/S"
anno = rbind(anno,c("G1","MFA1",1.44)) 
anno = rbind(anno,c("M/G1","DSE2",9.329))
anno = rbind(anno,c("M/G1","DSE1",7.869))
anno = rbind(anno,c("M/G1","CTS1",11.200))
anno$Type = factor(anno$Type,levels = c("M/G1","G1","G1/S","S","G2/M"))
anno = anno[order(anno$Type),]

anno$scer_gene = paste("SCER-",anno$Marker,sep="")


#factor(anno$Type,levels = c("M/G1","G1","S","G2/M")
#cell_cycle_df
cell_cycle_big_df = read.csv("data/cc//cell_cycle.csv")

cell_cycle_big_df = cell_cycle_big_df %>% inner_join(cell_cycle_df, by=c("ORF"="ORF"))
cell_cycle_big_df$NAME2  = paste("SCER-",cell_cycle_big_df$NAME,sep="")

#cell_cycle_df
#table(cell_cycle_big_df$Peak)

cell_cycle_big_df$Peak[cell_cycle_big_df$Peak == "S/G2"] = "G2/M"

cc_list = list(genes_name_trans=genes_name_trans,cell_cycle_big_df=cell_cycle_big_df,anno=anno)
#saveRDS(cc_list,"data/cc_list.RDS")
#saveRDS(gff_in, "data/gff_in.RDS")

haploid_genes = read.csv("data/haploid_genes.csv")
haploid_genes[!haploid_genes$name %in% genes_name_trans$gene_name,]
haploid_genes = haploid_genes %>% filter(name != "DPS2") 
#gene_name
MF2 = "MF(ALPHA)2"
MF1 = "MF(ALPHA)1"
idxs = which(!haploid_genes$name %in% genes_name_trans$gene_name)
haploid_genes$name[idxs][1] = MF1
haploid_genes$name[idxs][2] = MF2
haploid_genes$name[idxs][3] = "AFB1"
haploid_genes$name[idxs][4] = "MATALPHA2"
haploid_genes[which(!haploid_genes$name %in% genes_name_trans$gene_name),]
haploid_genes_df = haploid_genes %>% inner_join(genes_name_trans,by=c("name"="gene_name"))
scer_hsg_genes = paste0("SCER-",haploid_genes$name)
scer_hsg_genes_cell_cycle = scer_hsg_genes[scer_hsg_genes %in% cell_cycle_big_df$NAME2]
