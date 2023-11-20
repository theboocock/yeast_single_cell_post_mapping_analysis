library(Matrix)
# TODO: get working # 
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
for (folder in folders){
  cross_tmp = basename(folder)
  segdata = readRDS(glue("{folder}/segData.RDS"))
  umi=median(segdata$barcode.features$nUMI)
  n_expressed = median(colSums(ginformative_rds[[cross_tmp]]))
  n_cells = dim(segdata$Yr)[1]
  
  n_genes = dim(segdata$Yr)[2]
  
  if(cross_tmp %in% de_novo){
    tmp_mat = scale(t(segdata$Gsub))
    
    gs = ifelse(segdata$Gsub > .5,1,0)
    gs2 = hamming(gs)
    gs3 = 1-gs2/ncol(segdata$Gsub)
    
    
    diag(gs3) = 0
    idx = upper.tri(gs3)
    gs3[idx] = 0 
    idx_clones = which(gs3 > 0.8,arr.ind = T)
    non_unique = unique(c(idx_clones[,1],idx_clones[,2]))
    prop_unique = length(non_unique) / n_cells
    
    
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
  }
  out_df = rbind(out_df,data.frame(cross=cross_tmp, n_transcripts=n_genes, n_cells=n_cells, median_umi=umi, expressed_snps=n_expressed,prop_unique=prop_unique,folder=normalizePath(folder)) )
}

### ###
# 

#.parents)[c(8)],
#'16_GP3_2999_3000_3001_3049_May10'=names(crosses.to.parents)[c(9,11)],
#'17_GP4_3003_3043_3028_381_May10'=names(crosses.to.parents)[c(13,16)])
#good.dips$`18_3051_May10`

#table(unlist(good.dips))
#good.dips

out_df$type = "sceQTL"


base_dir = "../rproj/out/"
out_df2 = out_df

out_df = data.frame()
for(base_path in names(good.dips)){
  #print(folder) 
  print(base_path)
  in_dip= good.dips[[base_path]]
  ase_data = readRDS(paste0(base_dir,"/",base_path, "/aseData.RDS"))
  diploid_assignment = readRDS(paste0(base_dir,"/",base_path, "/diploid_assignments.RDS"))
  mat = ase_data$counts %>% as.dgCMatrix.spam()
  colnames(mat) = ase_data$barcodes  
  rownames(mat) = ase_data$features$V1
  for(dip in in_dip){
    print(dip)
    
    bc=diploid_assignment$barcode[diploid_assignment$diploid_name == dip]
    n_cells = length(bc)
    numi = median(colSums(mat[,bc]))
    nbin = readRDS(paste0(base_dir,"/",base_path,"/bbin_",dip,"_CCmanual.RDS"))
    ngenes=length(unique(nbin$gene))
    
    path = normalizePath(paste0(base_dir,"/",base_path))
    
    out_df = rbind(out_df, data.frame(cross=dip,n_transcripts=ngenes,n_cells=n_cells,median_umi=numi, expressed_snps=NA,prop_unique=NA,folder=path,type="ASE"))
    #dim(nbin)
    #mat[,bc]
    #base.dir = paste("../rproj/out/",base_path,dip,sep="")
    #list.files(base.dir)
    #readRDS(past)
  }
  #paste0(folder,"/",)
}

summary_table = rbind(out_df2,out_df)
#sum(summary_table$n_cells)
summary_table$type[c(5,8,9)] = "ASE_replication"
summary_table %>% filter(type == "ASE_replication")
#sum(summary_table$n_cells)
write_csv(summary_table,"Tables/s4.csv")