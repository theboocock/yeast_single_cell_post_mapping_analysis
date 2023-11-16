

### Dataset's used for analysis #### 

###
library(Matrix)
ginformative_rds = readRDS("../rproj/out/combined/ginformative.RDS")

folders = list.dirs("../rproj/out/combined/", full.names=T)
base.dir = "../rproj/out/combined/"
crosses = c("Ap","A","B","3004")
folders = paste(base.dir,crosses,sep="")

de_novo = c("A","B","3004")
#data.frame(crosses)

out_df = data.frame()


hamming <- function(X) {
  D <- (1 - X) %*% t(X)
  D + t(D)
}
library(glue)


unique_cells = data.frame(cross=c('A','B','3004'),n_unique_cells=c(1987,675,0))


for (folder in folders){
  cross_tmp = basename(folder)
  segdata = readRDS(glue("{folder}/segData.RDS"))
  umi=median(segdata$barcode.features$nUMI)
  n_expressed = median(colSums(ginformative_rds[[cross_tmp]]))
  n_cells = dim(segdata$Yr)[1]
  
  n_genes = dim(segdata$Yr)[2]
  
  if(cross_tmp %in% de_novo){
    n_unique_cells = unique_cells$n_unique_cells[unique_cells$cross == cross_tmp]
    #tmp_mat = scale(t(segdata$Gsub))
    prop_unique = n_unique_cells/n_cells
  #  gs = ifelse(segdata$Gsub > .5,1,0)
    #gs2 = hamming(gs)
    #gs3 = 1-gs2/ncol(segdata$Gsub)
    
    
    #diag(gs3) = 0
    #idx = upper.tri(gs3)
    #gs3[idx] = 0 
    #idx_clones = which(gs3 > 0.8,arr.ind = T)
    #non_unique = unique(c(idx_clones[,1],idx_clones[,2]))
    #prop_unique = length(non_unique) / n_cells
    
    
    #tmp_mat = t(tmp_mat)
    #tmp_mat = tcrossprod(tmp_mat)
    #tmp_mat = tmp_mat /(ncol(segdata$Gsub)-1)
    
    
    
    
    # diag(tmp_mat) = 0
    #idx  = upper.tri(tmp_mat)
    ## tmp_mat[idx] = 0 
    # idx_clones = which(tmp_mat > 0.7,arr.ind = T)
    # non_unique = unique(c(idx_clones[,1],idx_clones[,2]))
    # prop_unique = length(non_unique) / n_cells
  }else{
    prop_unique = NA
    n_unique_cells = NA
  }
  out_df = rbind(out_df,data.frame(cross=cross_tmp, n_transcripts=n_genes, n_cells=n_cells, median_umi=umi, expressed_snps=n_expressed,n_duplicated=n_unique_cells,prop_duplicated=prop_unique,folder=normalizePath(folder)) )
}


### ###
# 
#using 1-hamming distance of 0.8 as the threshold
#for BYxRM:
#  1987 non-unique
#27744 total
#0.072 fraction of non-unique
#for YPSxYJM
#675 non-unique
##44784 total
#0.015 fraction of non-unique
#.parents)[c(8)],
#'16_GP3_2999_3000_3001_3049_May10'=names(crosses.to.parents)[c(9,11)],
#'17_GP4_3003_3043_3028_381_May10'=names(crosses.to.parents)[c(13,16)])
#good.dips$`18_3051_May10`

#table(unlist(good.dips))
#good.dips

out_df$type = "eQTL"


base_dir = "../rproj/out/"
num_experiments = c(1,1,2,1)
num_libraries = c(2,2,4,2)
## update fraction no unique
#unique = c(NA,0.072,0.015,0)
out_df$n_experiments=num_experiments
out_df$n_libraries = num_libraries
out_df$cross_name = c("393 BYxRM","A","B","C")
out_df$parent1 = c("BY","BY","YJM145","CBS2888")
out_df$parent2 = c("RM","RM","YPS163","YJM981")
out_df$effect_directions = c("-BY,+RM","-BY,+RM","-YJM145,+YPS163","-CBS2888,+YJM981")
out_df2 = out_df

out_df = data.frame()
for(base_path in names(good.dips)){
  #print(folder) 
  print(base_path)
  in_dip= good.dips[[base_path]]
  ase_data = readRDS(paste0(base_dir,"/",base_path, "/aseData.RDS"))
  #phased_counts = readRDS(paste0(base_dir,"/","phasedCounts"))
  diploid_assignment = readRDS(paste0(base_dir,"/",base_path, "/diploid_assignments.RDS"))
  mat = ase_data$counts %>% as.dgCMatrix.spam()
  colnames(mat) = ase_data$barcodes  
  rownames(mat) = ase_data$features$V1
  for(dip in in_dip){
    phased_counts = readRDS(paste0(base_dir,"/",base_path,"/phasedCounts_",dip,".RDS"))
    #print(dip)
    
    bc=diploid_assignment$barcode[diploid_assignment$diploid_name == dip]
    n_cells = length(bc)
    numi = median(colSums(mat[,bc]))
    nbin = readRDS(paste0(base_dir,"/",base_path,"/bbin_",dip,"_CCmanual.RDS"))
    ngenes=length(unique(nbin$gene))
    #expressed_snps = median(colSums((phased_counts$par1.ASEcounts[bc,] + phased_counts$par2.ASEcounts[bc,]) > 0))
    
    expressed_snps = median(colSums((ase_data$rC[,diploid_assignment$diploid_name == dip] + ase_data$aC[,diploid_assignment$diploid_name == dip]) > 0))
    
    path = normalizePath(paste0(base_dir,"/",base_path))
    
    out_df = rbind(out_df, data.frame(cross=dip,n_transcripts=ngenes,n_cells=n_cells,median_umi=numi, expressed_snps=expressed_snps,n_duplicated=NA,prop_duplicated=NA,folder=path,type="ASE"))
    #dim(nbin)
    #mat[,bc]
    #base.dir = paste("../rproj/out/",base_path,dip,sep="")
    #list.files(base.dir)
    #readRDS(past)
  }
  #paste0(folder,"/",)
}
out_df$n_experiments = 1 
out_df$n_libraries = 1
out_df$cross_name = c("A","B","C","B","C","A")
out_df$parent1 = c("BY","YJM145","CBS2888","YJM145","CBS2888","BY")
out_df$parent2 = c("RM","YPS163","YJM981","YPS163","YJM981","RM")
out_df$effect_directions = c("-BY,+RM","-YJM145,+YPS163","-CBS2888,+YJM981","-YJM145,+YPS163","-CBS2888,+YJM981","-BY,+RM")
out_df$type = c("ASE_REP","ASE","ASE","ASE_REP","ASE_REP","ASE")
#out_df
#out_df3 = out_df[c(2,3,6),]
#out_df3

#out_df$cross_ame = 

#out_df$

summary_table = rbind(out_df2,out_df)
#sum(summary_table$n_cells)
#summary_table$type[c(5,8,9)] = "ASE_replication"
#summary_table %>% filter(type == "ASE_replication")
#sum(summary_table$n_cells)
write_tsv(summary_table,"tables//s4.csv")
openxlsx::write.xlsx(summary_table,"tables/s4.xlsx")
