
## Load in data for TopGo
gff_in= rtracklayer::import.gff("/media/theboocock/Data/Dropbox/PHDTHESIS/Single-CellRNASEQ/eqtls/ref_data/yeast/gffs/saccharomyces_cerevisiae.gff")
gff_in = gff_in %>% filter(type == "gene")
gene_ontologies = apply(as.data.frame(gff_in),1, function(x){x$Ontology_term})
names(gene_ontologies) = gff_in$Name
### Filter for genes used in the anlysis #
GO2geneID = inverseList(gene_ontologies)

##
gff_in2 = (gff_in %>% filter(orf_classification != "Dubious"))
genes_all = gff_in2$Name
get_enrichment = function(gene_list, background=NULL){
  #print(x)
  #if(!is.null(background)){
  #  background = NULL
  #}
  gene_ontologies_tmp = (gene_ontologies[names(gene_ontologies) %in% background])
  all_genes = as.factor(rep(0,length(gene_ontologies_tmp)))
  genes_all_tmp = gene_ontologies[names(gene_ontologies) %in% background]
  g1 = gene_ontologies_tmp[names(gene_ontologies_tmp) %in% background]
  g1_f = factor(as.numeric(names(g1) %in% gene_list))
  names(g1_f) = names(g1)
  GOData = new("topGOdata",ontology="BP",allGenes=g1_f,annot=annFUN.gene2GO,gene2GO = gene_ontologies)
  test.stat =new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFis <- runTest(GOData, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOData, classic = resultFis, orderBy = "weight", ranksOf = "classic", topNodes = 300)
  
  bp_res = allRes 
  bp_res$log2FC = log2(allRes$Significant / allRes$Expected + 1)
  #annot=annFUN.gene2GO, gene2GO = gene_ontologies)
  GOData = new("topGOdata",ontology = "MF", allGenes=g1_f, annot=annFUN.gene2GO, gene2GO = gene_ontologies)
  test.stat =new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFis <- runTest(GOData, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOData, classic = resultFis, orderBy = "weight", ranksOf = "classic", topNodes = 300)
  mf_res = allRes
  mf_res$log2FC = log2(allRes$Significant / allRes$Expected + 1)
  
  mf_res$type = "MF"
  bp_res$type = "BP"
  gg  = rbind(mf_res, bp_res)
  #print(gg)
  gg$classic = str_replace_all(gg$classic,"<1e-30","1e-30")
  gg$classic_adj = p.adjust(gg$classic,method = "BH")
  out_mf_bp = rbind(mf_res, bp_res)  %>% mutate(classic = as.numeric(str_replace_all(classic,"<",""))) %>% filter(log2FC > 2 & classic < 1e-5)
  
  #keytypes(org.Sc.sgd.db)
  
  aaa = bitr(names(g1_f),fromType = "ORF",toType = "ENTREZID",OrgDb = org.Sc.sgd.db)
  
  keg = enrichKEGG(names(g1_f)[g1_f == 1],universe = aaa$ORF,organism="sce")
  keg_df = keg@result %>% filter(qvalue < 0.05)
  return(list(out_mf_bp=out_mf_bp,kegg=keg_df,out_mf_bp_raw=gg))
}
#get_rnaseq_enrichment = function(gene_list, background=NULL){

#  )  
#}

whdquantile <- function(x, probs, weights = NA) {
  cdf.gen <- function(n, p) return(function(cdf.probs) {
    pbeta(cdf.probs, (n + 1) * p, (n + 1) * (1 - p))
  })
  wquantile.generic(x, probs, cdf.gen, weights)
}

wquantile.generic <- function(x, probs, cdf.gen, weights = NA) {
  n <- length(x)
  if (any(is.na(weights)))
    weights <- rep(1 / n, n)
  nw <- sum(weights)^2 / sum(weights^2) # Kish's effective sample size
  
  indexes <- order(x)
  x <- x[indexes]
  weights <- weights[indexes]
  
  weights <- weights / sum(weights)
  cdf.probs <- cumsum(c(0, weights))
  
  sapply(probs, function(p) {
    cdf <- cdf.gen(nw, p)
    q <- cdf(cdf.probs)
    w <- tail(q, -1) - head(q, -1)
    sum(w * x)
  })
}


hotspot_enrichment_and_function = function(cross,hotspot_peaks,combined_peaks,lod_df_beta, cell_cycle_beta_df2,segdata,fdr_fx,background,out_name=NULL){
  #print(folder)
  #cross = basename(folder)
  peaks = hotspot_peaks %>% dplyr::filter(cell_cycle_stage == "combined") %>% filter(in.hotspot)  %>% dplyr::group_by(bin) %>% dplyr::summarise(n=dplyr::n()) %>% arrange(-n)
  sig_hotspot = qpois(1-(.05/nrow(peaks)),ceiling(mean(peaks$n))) + 1
  
  
  parents =crosses.to.parents[[cross]]
  causal_hit= causal_hits%>% filter(shared.parent  %in% parents)  %>% mutate(seqnames=chr) %>% as_granges()
  #causal_hit = causal_hit_3004 %>% mutate(seqnames=chr) %>% as_granges()
  qtl = qtl_bloom %>% filter(cross == !!cross)
  vcf_names = str_remove_all(parents,pattern="a|x")
  vcf_annotation = vcf_annotation_16 %>% filter((grepl(vcf_names[1],ref_strains) & grepl(vcf_names[2],alt_strains)) |  (grepl(vcf_names[2],ref_strains) & grepl(vcf_names[1],alt_strains)))
  ## Missense of STOP
  size = unlist(lapply(str_split(vcf_annotation$mutation,"\\|"), function(x){x[3]}))
  vcf_annotation$impact  = size 
  #vcf_annotation = vcf_annotation[size %in% c("HIGH","MODERATE")]
  # Get significant peaks ##
  peak_names_sig =peaks 
  str_split_peaks = str_split(peak_names_sig$bin,':')
  chrom = unlist(lapply(str_split_peaks, function(x){x[1]}))
  left = unlist(lapply(str_split_peaks, function(x){x[2]}))
  start = as.numeric(unlist(lapply(str_split(left,"-"), function(x){x[1]})))
  end = as.numeric(unlist(lapply(str_split(left,"-"), function(x){x[2]})))
  dfr = data.frame(seqnames=chrom,start=start, end=end, bin=peak_names_sig$bin,n=peak_names_sig$n)  %>% as_granges()
  new_bins  = data.frame(seqnames=chrom,start=start, end=end, bin=peak_names_sig$bin,n=peak_names_sig$n) %>% as_granges() %>% reduce_ranges() %>% join_overlap_inner(dfr) %>% as_data_frame()
  
  cc_lods = lod_df_beta
  #hotspot_peaks
  
  #right = lapply(str_split_peaks,function(x){x[3]})
  
  
  #cc_lods = lod_dfs[[cross]]
  #cc_red = cc_lods %>% filter(seqnames != "chrIII") %>% reduce_ranges() 
  #cc_red_df = cc_red %>% as_data_frame() %>% mutate(bin=paste(seqnames,paste(start,end,sep="-"),sep=":")) %>% as_granges()
  #cc_red$bin=cc_red_df$bin
  #icc_red$snp = cc_red_df
  #cc_red_big = cc_red  %>% join_overlap_inner(cc_lods)
  #for(region in unique(cc_red$bin)){
  #  print(region)
  ##  cc_red_tmp = cc_red  %>% filter(bin == region)
  #  snpid = cc_red_tmp$
  
  #}
  #lod_df_beta
  
  ##peak_names_sig$bin
  hotspot_peaks_tmp = hotspot_peaks %>% filter(bin %in% peak_names_sig$bin)
  for(i in 1:nrow(new_bins)){
    old_bin = new_bins$bin[i]
    new_bin = paste(new_bins$seqnames[i],paste(new_bins$start[i],new_bins$end[i],sep="-"),sep=":")
    hotspot_peaks_tmp = hotspot_peaks_tmp %>% mutate(bin = ifelse(bin == old_bin,new_bin,bin))
  }
  hotspot_str_list = c()
  annotation_list = list()
  j= 1    
  betas_pheno = load_pheno_maps(cross, segdata)
  traits = rownames(betas_pheno)
  # TODO test tommorrow
  str_split_2 = str_split(unlist(lapply(str_split(hotspot_peaks_tmp$bin,":"), function(x){x[2]})),"-")
  hotspot_pos = round((as.numeric(unlist(lapply(str_split_2,function(x){x[1]}))) +  as.numeric(unlist(lapply(str_split_2,function(x){x[2]}))))/2)
  hotspot_peaks_tmp$hotspot_pos = hotspot_pos
  #hotspot_out_all_confidence_interval = data.frame()
  for(hotspot in unique(hotspot_peaks_tmp$bin)){
    #print(hotspot)
    tmp_df = hotspot_peaks_tmp[hotspot_peaks_tmp$bin == hotspot,]
    #quants = quantile(tmp_df$pos,c(0.05,.1,.2,.8,.9,.95))
    #print(quants)
    
    #segdata
    buffer = 5e3
    g  = hotspot_peaks_tmp %>%  filter(bin ==hotspot) %>% arrange(-LOD) #filter(LOD > 10) %>% arrange(-LOD) 
    #g %>% filter(LOD > 10)
    g2 = g %>% arrange(-LOD) %>% head(n=20)
    #print(g2)
    pos_one = quantile(g2$CI.l,prob=.1,F,T,3)
    pos_two = quantile(g2$CI.r,prob=0.9,F,T,3)
    #whdquantile(g$CI.l,prob=.1,weights=g$LOD)
    #whdquantile(g$CI.r,prob=.9,weights=g$LOD)
    
    #print(g)
    chrom = unlist(lapply(str_split(colnames(segdata$Gsub),"_"), function(x){x[1]}))
    pos = as.numeric(unlist(lapply(str_split(colnames(segdata$Gsub),"_"), function(x){x[2]})))
    #print(chrom)
    ### Yea of course that doesn't work because I have to find the nearest marker ###
    chrom_tmp = chrom[chrom == g$chrom[1]]
    pos_tmp = pos[chrom == g$chrom[1]]
    
    
    idx_one = which(min(abs(pos_tmp - pos_one)) == abs(pos_tmp - pos_one))
    idx_two = which(min(abs(pos_tmp - pos_two)) == abs(pos_tmp - pos_two))
    #print(idx_one)
    #print(idx_two)
    pos_one = pos_tmp[idx_one]
    pos_two = pos_tmp[idx_two] 
    left = which(g$chrom[1]==chrom & pos == pos_one) 
    right = which(g$chrom[1]==chrom & pos == pos_two)# +2
    #print(left) 
    #print(right)
    
    #print(idx_one)
    #print(idx_two)
    ##left = idx_one - 2
    #right = idx_two + 2
    #left  = which(chrom == g$chrom[1] & pos == g$CI.l[1]) -2
    #right =  which(chrom == g$chrom[1] & pos == g$CI.l[1]) + 2
    
    chrom_l = chrom[left]
    chrom_r = chrom[right]
    pos_l = pos[left]
    pos_r = pos[right]
    #print(left)
    if(chrom_l != g$chrom[1]){
      chrom_l = g$chrom[1]
      pos_l = 1
    }
    if(chrom_r != g$chrom[1]){
      chrom_r = g$chrom[1]
      pos_r = chrom_lengths$lengths[chrom_lengths$chrom == g$chrom[1]]
    }
    
    start= pos_l
    end = pos_r
    region_df = data.frame(seqnames=unique(g$chrom),start=start,end=end) %>% as_granges()
    ##
    #print(region_df)    
    genes_in_region = join_overlap_inner(gff_in2,region_df)
    
    
    #end(genes_in_region)
    genes_in_region$tpos = round((start(genes_in_region) +    end(genes_in_region))/2)
    #split_snp_id = str_split(colnames(all_pheno_maps[[cross]]),"_")
    #chrom_snps  =lapply(split_snp_id, function(x){x[1]})
    #pos_snps = as.numeric(lapply(split_snp_id, function(x){x[2]}))
    #snp_ids_in_pheno = paste(chrom_snps,pos_snps,sep="_")
    #snp_ids2 = paste(g$chrom[1], g$pos[1],sep="_")
    #cc_lods = lod_dfs[[cross]]
    dd = data.frame()
    
    out_cc = cell_cycle_beta_df2[which(cell_cycle_beta_df2$marker == g$peak.marker[1]),]
    
    # (beta_lists[["3004"]]$beta_list$G1[which(snp_ids_in_pheno %in% snp_ids2)])
    # which(snp_ids_in_pheno %in% snp_ids2)
    # = all_pheno_maps[[cross]][, which(snp_ids_in_pheno %in% snp_ids2)]
    betas_pheno_tmp = betas_pheno[,colnames(betas_pheno) %in% g$peak.marker[1]]
    lods = get_lod_from_r(betas_pheno_tmp, n=dim(seg.recoded[[cross]])[2])
    df_qtl = data.frame(traits=traits,betas=betas_pheno_tmp,lods=lods)
    #betaswhich(snp_ids_in_pheno %in% snp_ids2)
    #hotspot_id 
    gene_list = g$transcript
    background = colnames(segdata$Yr)
    background = (background[background %in% gff_in$Name[as.character(seqnames(gff_in)) != g$chrom[1]]])
    
    enrich_df = get_enrichment(gene_list = g$transcript,background = background)
    #lod_dfs$`out/combined//3004`
    cell_cycle_causals = region_df %>% join_overlap_inner(cc_lods) %>% as_data_frame()
    causals_round_robin = region_df%>% join_overlap_inner(causal_hit) %>% as_data_frame()
    assoc_round_robin = region_df %>% join_overlap_inner(qtl) %>% as_data_frame()
    annotation_provean_snps = region_df %>% join_overlap_inner(vcf_annotation) %>% as_data_frame() # %>% write_tsv("candidates.tsv")
    no_transcripts = nrow(g)
    #out_file = glue("out/fine_mapping/{cross}_{hotspot}_{no_transcripts}_finemapping.xlsx")
    if(!is.null(out_name)){
      dir.create(glue("out/fine_mapping/{out_name}"),recursive = T)
      out_file = glue("out/fine_mapping/{out_name}/{hotspot}_{no_transcripts}_finemapping.xlsx")
      #print(out_file)
    }else{
      dir.create(glue("out/fine_mapping/{cross}"),recursive = T)
      out_file = glue("out/fine_mapping/{cross}/{hotspot}_{no_transcripts}_finemapping.xlsx")
      #print(out_file)
    }
    hotspot_str = glue("{hotspot}_{no_transcripts}")
    if(nrow(cell_cycle_causals) == 0){
      #print(nrow(cell_cycle_causals))
      cell_cycle_causals= data.frame(DATA="NONE")
    }
    if(nrow(causals_round_robin) == 0){
      #print(nrow((causals_round_robin)))
      causals_round_robin = data.frame(DATA="NONE")
    }
    if(ncol(df_qtl) == 0){
      df_qtl = data.frame(DATA="NONE")
    }
    
    if(nrow(out_cc) == 0){
      out_cc = data.frame(DATA="NONE")
    }
    ab = genes_in_region %>% as_data_frame() %>% dplyr::select(-Note,-Ontology_term,-Parent,-Alias) 
    if(nrow(ab) == 0){
      ab = data.frame(DATA="NONE")
    }
    #ab %>% write.table("test.txt")
    
    peak_merged = hotspot_peaks[hotspot_peaks$transcript %in% g$transcript &
                                            hotspot_peaks$peak.marker %in% g$peak.marker,]
    
    peak_merged = peak_merged %>% arrange(-LOD)
    peak_merged = peak_merged %>% inner_join(genes_name_trans,by=c("transcript"="gene_id"))
    if(!is.null(fdr_fx)){
    local_eqtls = combined_peaks %>% mutate(cis_fdr=fdr_fx$c.rtoFDRfx(LOD)) %>% filter(cis_fdr < 0.05)%>% filter(chrom == g$chr[1] & tchr == g$chr[1])
    }else{
      local_eqtls = combined_peaks %>% filter(chrom == g$chr[1] & tchr == g$tchr[1])
      
    }
    # TODO CI doesn't need to be within the Trans CI that's weird. fix it
    
    local_eqtls = local_eqtls %>% filter(tpos > start & tpos < end)
    local_eqtls = local_eqtls %>% inner_join(genes_name_trans,by=c("transcript"="gene_id"))
    local_eqtls = local_eqtls %>% arrange(-LOD)
    if(nrow(local_eqtls) ==0){
      local_eqtls = data.frame(DATA="NONE")
    }
    hotspot_str_list = c(hotspot_str_list, hotspot_str)
    annotation_list[[j]]=list(cell_cycle_causals=cell_cycle_causals,cell_cycle_assoc=dd,causals_round_robin=causals_round_robin,
                                    assoc_round_robin=assoc_round_robin, complete_assoc=df_qtl, provean_snps = as.data.frame(annotation_provean_snps),
                                    GO=enrich_df$out_mf_bp, KEGG=enrich_df$kegg,GO_raw=enrich_df$out_mf_bp_raw,directional_hotspot_distal=peak_merged,
                                   directional_hotspot_local=local_eqtls,genes_in_region=ab,eqtl_region=region_df)
    
    #annotation_list[[j]]
    
    #region_df = region_df %>% as_data_frame()
    #content = openxlsx::createWorkbook(out_file)
    #openxlsx::
    #list()
    
    openxlsx::write.xlsx(annotation_list[[j]],file=out_file)
    
   # hotspot_out_all_confidence_interval = rbind(hotspot_out_all_confidence_interval,region_df)
    #openxlsx::write.xlsx(cell_cycle_causals, file=out_file, sheetName = "cell_cycle_causals")
    #openxlsx::write.xlsx(out_cc,file=out_file, sheetName = "cell_cycle_assoc",overwrite =F)
    #openxlsx::write.xlsx(causals_round_robin,file=out_file,sheetName = "causals_round_robin",overwrite = F)
    #openxlsx::write.xlsx(assoc_round_robin,file=out_file,sheetName = "assoc_round_robin",overwrite = F)
    #openxlsx::write.xlsx(df_qtl, file=out_file, sheetName = "complete_assoc", overwrite = F)
    #openxlsx::write.xlsx(as.data.frame(annotation_provean_snps),file=out_file, sheetName = "provean_snps",overwrite = F)
    #openxlsx::write.xlsx(enrich_df$out_mf_bp, file=out_file, sheetName="GO",overwrite = F)
    #openxlsx::write.xlsx(enrich_df$kegg, file=out_file, sheetName="KEGG", overwrite = F)
    #openxlsx::write.xlsx(peak_merged,file=out_file,sheetName="directional_hotspot_distal",overwrite = F)
    #openxlsx::write.xlsx(local_eqtls, file=out_file,sheetName="directional_hotspot_local",overwrite = F)
    #openxlsx::write.xlsx(ab, file=out_file, sheetName = "genes_in_hotspot",overwrite = F)
    #openxlsx::write.xlsx(region_df, file=out_file,sheetName = "region_df", overwrite = F)
    j =j +1
  }
  names(annotation_list) = hotspot_str_list 
  hotspots_list=hotspot_peaks_tmp
  #i=i+1
  return(list(annotation_list=annotation_list,hotspot_list=hotspots_list,hotspot_threshold=sig_hotspot))
}


