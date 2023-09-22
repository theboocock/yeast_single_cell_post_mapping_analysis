get_theta=function(beta,se){
  
  weights = 1/se^2
  top = sum(beta * weights)
  bot = sum(weights)
  
  return(top/bot)
  
}
get_cochran_q = function(theta,beta,se){
  weights = 1/se^2
  # print(weights)
  Q = sum(weights* (theta-beta)^2)
  return(Q) 
}

get_lod_from_r = function(r,n=365){
  return(-n * (log(1- r^2) / (2 *log(10))))
}

make_cc_lod_traces = function(cell_cycle_df){
  
  cell_cycle_rds_df %>% ggplot(aes(y=lod,x=pos)) + geom_point() + facet_grid( rows = vars(cell_cycle),  cols=vars(chrom),scales="free",) +
    theme_classic() + theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=1))
}

convert_chrom_to_simple_factor = function(chrom,reverse=F){
  new_chrom = str_replace_all(chrom,"chr","")
  order_levels = c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")
  #rev_order_levels = rev(order_levels)
  if(!reverse){
    new_chrom = factor(new_chrom, levels=order_levels)
  }else{
    new_chrom = factor(new_chrom, levels=rev(order_levels))
  }
  #hotspot_peaks$chrom_short_f = factor(hotspot_peaks$chrom_short,levels=order_levels)
  return(new_chrom)
}


convert_cell_cycle_to_simple_factor = function(cell_cycle){
  cell_cycle = as.character(cell_cycle)
  order_levels = c("M/G1","G1","G1/S","S","G2/M")
  if(any(str_count(cell_cycle,":") > 0 )){
    cell_cycle = str_replace_all(cell_cycle,":","/")
  }
  cell_cycle = factor(cell_cycle,levels=order_levels)
  return(cell_cycle)
}


load_pheno_maps = function(cross,segdata){
  #all_pheno_maps = list()
  #for(folder in folders){
  pheno_and_genotype = list(seg_recoded=seg.recoded, pheno_average = pheno_extracted)
  seg_tmp =t(pheno_and_genotype$seg_recoded[[cross]])
    # seg tmp
  seg_reduced = seg_tmp[,which(colnames(seg_tmp) %in% colnames(segdata$Gsub))]
    #plot(colMeans(seg_reduced))
  seg_reduced = apply(seg_reduced, 2, scale)
  pheno_resid = pheno_and_genotype$pheno_average[[cross]]
    # pheno resid 
  scaled_pheno = apply(pheno_resid,2, scale)
  all_pheno = crossprod(scaled_pheno,seg_reduced)/(nrow(seg_reduced) - 1)
  return(all_pheno)   
}

make_cc_lod_df = function(cell_cycle_rds, experiment){
  df= cell_cycle_rds %>% as.data.frame() %>%
    rownames_to_column(var = "cell_cycle") %>%  
    filter(!grepl("seurat",cell_cycle)) %>%
    pivot_longer(!cell_cycle,names_to="marker",values_to="lod")
  df$cell_cycle = str_replace_all(str_replace_all(df$cell_cycle,"cell_cycle",""),"/",":")
  split_str = str_split(df$marker,"_")
  chrom = unlist(lapply(split_str, function(x){x[1]}))
  pos = as.numeric(unlist(lapply(split_str, function(x){x[2]})))
  ref = (unlist(lapply(split_str, function(x){x[3]})))
  alt = unlist(lapply(split_str, function(x){x[4]}))
  df$chrom = chrom
  df$pos = pos
  df$ref = ref
  df$alt = alt
  df$cell_cycle = factor(df$cell_cycle)
  df$chrom = factor(df$chrom)
  df$experiment = as.character(experiment)
  return(df)
  #return(df)
}

make_seurat_lod_df = function(cell_cycle_rds, experiment){
  df= cell_cycle_rds %>% as.data.frame() %>%
    rownames_to_column(var = "cell_cycle") %>%  
    filter(grepl("seurat",cell_cycle)) %>%
    pivot_longer(!cell_cycle,names_to="marker",values_to="lod")
  #df$cell_cycle = str_replace_all(str_replace_all(df$cell_cycle,"seura",""),"/",":")
  split_str = str_split(df$marker,"_")
  chrom = unlist(lapply(split_str, function(x){x[1]}))
  pos = as.numeric(unlist(lapply(split_str, function(x){x[2]})))
  ref = (unlist(lapply(split_str, function(x){x[3]})))
  alt = unlist(lapply(split_str, function(x){x[4]}))
  df$chrom = chrom
  df$pos = pos
  df$ref = ref
  df$alt = alt
  df$cell_cycle = factor(df$cell_cycle)
  df$chrom = factor(df$chrom)
  df$experiment = as.character(experiment)
  return(df)
  #return(df)
}

getBinAssignment=function(cPf, cbin50k,markerGR) {
  tlinks=makeGRangesFromDataFrame(data.frame(chr=cPf$chr, start=cPf$pos, end=cPf$pos, strand="*")) 
  
  tlink.cnt=findOverlaps(tlinks, unlist(cbin50k))
  idx_no_match = which(is.na((match(1:length(tlinks),queryHits(tlink.cnt) ))))
  cmax=sapply(sapply(split(markerGR, seqnames(markerGR)), start), max)
  cPf = correctPeakPositions(cPf,idx_no_match,cmax)
  tlinks=makeGRangesFromDataFrame(data.frame(chr=cPf$chr, start=cPf$pos, end=cPf$pos, strand="*")) 
  
  tlink.cnt=findOverlaps(tlinks, unlist(cbin50k))
  ub=unlist(cbin50k)
  ub.id=paste0(seqnames(ub), ':', start(ub), '-', start(ub)+width(ub)-1)
  cPf$bin=ub.id[subjectHits(tlink.cnt)]
  return(cPf)
}

correctPeakPositions = function(cPf, idx_no_match,cmax){
  for(idx in idx_no_match){
    cPf$pos[idx]=cmax[names(cmax) == cPf$chrom[idx]]
  }
  return(cPf)
}

get_granges_hotspot = function(x){
  chrom = unlist(lapply(str_split(names(table(x$bin)), ":"), function(x){x[1]}))
  left = as.numeric(unlist(lapply(str_split(lapply(str_split(names(table(x$bin)), ":"), function(x){x[2]}),"-"), function(x){x[1]})))
  right = as.numeric(unlist(lapply(str_split(lapply(str_split(names(table(x$bin)), ":"), function(x){x[2]}),"-"), function(x){x[2]})))
  granges_df = data.frame(seqnames=chrom,start=left,end=right) %>% as_granges()
  granges_df$bin = names(table(x$bin))
  return(granges_df)
}

makeBinTable_hotspot_cc = function(x,cbin50k){
  tlinks=makeGRangesFromDataFrame(data.frame(chr=x$chrom, start=x$pos, end=x$pos,strand="*")) 
  tlinks$interaction = x$has_interaction 
  tlinks_int = tlinks[tlinks$interaction]
  tlinks_no_int = tlinks[!tlinks$interaction]
  tlink.cnt=countOverlaps(unlist(cbin50k), tlinks_int)
  bin.table.int=data.frame(chr=seqnames(unlist(cbin50k)),
                       pos=round(start(unlist(cbin50k))+width(unlist(cbin50k))/2),
                       count=tlink.cnt)
  
  #tlinks_int = tlinks[tlinks$interaction]
  #tlinks_no_int = tlinks[!tlinks$interaction]
  tlink.cnt=countOverlaps(unlist(cbin50k), tlinks_no_int)
  bin.table.noint=data.frame(chr=seqnames(unlist(cbin50k)),
                           pos=round(start(unlist(cbin50k))+width(unlist(cbin50k))/2),
                           count=tlink.cnt)
  
  bin.table.int$no_int_count = bin.table.noint$count
  bin.table.int$total_count = bin.table.int$count + bin.table.int$no_int_count
  bin.table.int$int_count = bin.table.int$count
  bin.table.int$count = NULL
  return(bin.table.int)
  
}

makeBinTable_hotspot = function(x,cbin50k){
  tlinks=makeGRangesFromDataFrame(data.frame(chr=x$chrom, start=x$pos, end=x$pos, strand="*")) 
  tlink.cnt=countOverlaps(unlist(cbin50k), tlinks)
  bin.table=data.frame(chr=seqnames(unlist(cbin50k)),
                       pos=round(start(unlist(cbin50k))+width(unlist(cbin50k))/2),
                       count=tlink.cnt)
  return(bin.table)
  
}
getMarkerGRanges=function(g.counts) {
  
  g2m=tstrsplit(rownames(g.counts[[1]]),'_')
  markerGR=makeGRangesFromDataFrame(data.frame(chr=g2m[[1]], 
                                               start=as.numeric(g2m[[2]]),
                                               end=as.numeric(g2m[[2]]), 
                                               strand="*"))
  markerGR$name=rownames(g.counts[[1]])
  return(markerGR)
}


makeBinTable=function(cPf, cbin50k) {
  tlinks=makeGRangesFromDataFrame(data.frame(chr=cPf$chr, start=cPf$pos, end=cPf$pos, strand="*")) 
 
  #cPf[idx_no_match,]
  tlink.cnt=countOverlaps(unlist(cbin50k), tlinks)
  bin.table=data.frame(chr=seqnames(unlist(cbin50k)),
                       pos=round(start(unlist(cbin50k))+width(unlist(cbin50k))/2),
                       count=tlink.cnt)
  return(bin.table)
}

getcbin50k_genome = function(){
  
  cmin = rep(0,16)
  cmax = chrom_lengths$lengths[1:16]  
  names(cmin) = chrom_lengths$chrom[1:16]
  names(cmax) = chrom_lengths$chrom[1:16]
  cbin = data.frame(chr=chrom_lengths$chrom[1:16],start=cmin,end=cmax,strand="*")
  cbin = makeGRangesFromDataFrame(cbin)
  cbin50k = GenomicRanges::tile(cbin,width=50000)
  return(cbin50k)
}

getcbin50k= function(Gsub){
  markerGR=getMarkerGRanges(list(t(Gsub)))
  
  
  cmin=sapply(sapply(split(markerGR, seqnames(markerGR)), start), min)
  cmin[!is.na(cmin)]=0
  cmax=sapply(sapply(split(markerGR, seqnames(markerGR)), start), max)
  
  cbin=data.frame(chr=names(cmin),
                  start=cmin,
                  end=cmax, 
                  strand="*")
  cbin=makeGRangesFromDataFrame(cbin)
  ###########################################################################################################3
  # analyze within each sub experiment 
  cbin50k=GenomicRanges::tile(cbin, width=50000)
  return(cbin50k)
}

find_qtl = function(cc_rds_group,lod_threshold=5,lod_drop=2){
  idx_m = which(max(cc_rds_group$lod) == cc_rds_group$lod)[1]
  tmp_row = cc_rds_group[idx_m,]
  if(tmp_row$lod < lod_threshold){
    return(NULL)
  }
  lod = tmp_row$lod
  # to the right
  # to the left
  lod_drop_thresh = lod - lod_drop
  #print(cc_rds_group$lod[idx_m:1])
  left = which(cc_rds_group$lod[idx_m:1] < lod_drop_thresh)[1]  -1 
  right = which(cc_rds_group$lod[idx_m:length(cc_rds_group$lod)] < lod_drop_thresh)[1] -1
  #print(idx_m - left)
  #print(right + idx_m)
  #print(tmp_row)
  #print(cc_rds_group[(idx_m - left):(right+idx_m),])
  #print(cc_rds_group[idx_m - left,])
  #print(cc_rds_group[right + idx_m,])
  tmp_row$ci_left = cc_rds_group$pos[idx_m - left]
  tmp_row$ci_right = cc_rds_group$pos[right + idx_m]
  #if(tmp_row$ci_left){
  #  
  #}
  #if(tmp_row$ci_right){
  
  #}
  #print(tmp_row)
  return(tmp_row)
}

lod_drop_cc = function(cell_cycle_rds_df,lod_threshold=5,lod_drop=2){
  
  group_df = cell_cycle_rds_df %>% group_by(chrom, cell_cycle)
  gr_keys = group_keys(group_df)
  ## per cc map ##
  gr_map = group_df %>% group_map(~find_qtl(.x,lod_threshold =lod_threshold, lod_drop=lod_drop))
  ## total cc map ##
  all_cc = cell_cycle_rds_df %>% group_by(marker) %>% summarise(lod=sum(lod),chrom=chrom[1],pos=pos[1]) 
  keys_str= paste(gr_keys$chrom,gr_keys$cell_cycle, sep="_")
  names(gr_map) = keys_str
  out_df = bind_rows(gr_map,.id = "comparison")
  out_df$cell_cycle = unlist(lapply(str_split(out_df$comparison,"_"), function(x){x[2]}))
  out_df$chrom = unlist(lapply(str_split(out_df$comparison,"_"), function(x){x[1]}))
  
  return(out_df)
  # `group_all = all_cc %>% group_by(chrom) 
  #  all_map = group_all %>% group_map(~find_qtl(.x,lod_threshold = 5,lod_drop = 2))
  # names(all_map) = keys_str
  #bind_rows(all_map,.id = "comparison")
}
