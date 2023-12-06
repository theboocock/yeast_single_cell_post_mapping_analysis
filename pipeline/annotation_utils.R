# Load annotaiton matrix 


#in_dir = "/media/theboocock/Data/Dropbox/PHDTHESIS/projects/single_cell_2021/data/"

causal_hits = read.csv(glue("data/bloom_causal_variant_mapping.csv"), header=T)
#causal_hits = causal_hits %>% filter(shared.parent %in% c("YJM981x","CBS2888a"))
qtl_bloom = read.csv(glue("data/qtl_bloom.csv"))
pos_l = as.numeric(unlist(lapply(str_split(qtl_bloom$X1.5.LOD.drop.CI..left,"_"), function(x){x[2]})))
pos_r = as.numeric(unlist(lapply(str_split(qtl_bloom$X1.5.LOD.drop.CI..right,"_"), function(x){x[2]})))
qtl_bloom$pos_l = pos_l
qtl_bloom$pos_r =  pos_r
qtl_bloom = qtl_bloom %>% mutate(seqnames=chr,start=pos_l,end=pos_r) %>% as_granges()
vcf_annotation_16 = read_tsv("data/vcf/all_no_filt.txt",col_names=F) 
colnames(vcf_annotation_16) = c("chrom","pos","ref","alt","ref_strains","missing_strains","alt_strains","af","mutation")
split_annotation = str_split(vcf_annotation_16$mutation,"\\|")
#change  = lapply(split_annotation, function(x){x[11]})
#change = unlist(lapply(str_split(change, "\\."), function(x){x[2]}))
#original = substr(change, 1,3)
#pos = str_extract(change, "[0-9]+")
#aa_split = str_split(change,"[0-9]+")
#original = unlist(lapply(aa_split, function(x){x[1]}))
#change = unlist(lapply(aa_split, function(x){x[2]}))

#a(change)

#b  = "YKL198C"

#gene= unlist(lapply(split_annotation, function(x){x[5]}))

#provean_sc = c()

#for(i in 1:length(split_annotation)){
#  gene_tmp = gene[i]
#  change_tmp = change[i]
#  pos_tmp = pos[i]
#  change_tmp2 = a(change_tmp)#

# pr = provean_scores[[gene_tmp]][as.numeric(pos_tmp),change_tmp2]
#  provean_sc = c(provean_sc,pr)
#}

#02
#provean.score.file=paste0(working.dir, '../data/provean/') 

#load(glue("{in_dir}/provean/provean_scores_new.RData"))
#load("../data/provean/provean_scores_new.RData")
#meltedProveanScores=lapply(provean_scores, 
#                           function(x) {
#                             y=reshape::melt(x)
#                             names(y)[c(1,2)]=c('REF', 'ALT', 'score')
#                             y$POS=rep(1:nrow(x), ncol(x))
#                             return(y)
#                           })
#meltedProveanScores=data.table::rbindlist(meltedProveanScores, idcol='gene')
#meltedProveanScores$lookup=paste0(meltedProveanScores$gene, ':', meltedProveanScores$POS, ':', meltedProveanScores$ALT)


vcf_in = readRDS(glue("data/provean/1002codingEffects.RDS"))
#vcf_in = readRDS("../data/provean/1002codingEffects.RDS")

vcf_in$chrom = as.vector(seqnames(vcf_in$X))
vcf_in$start = start(ranges(vcf_in$X))




alts = str_split(str_remove_all(str_remove_all(str_remove_all(vcf_annotation_16$alt,pattern = "\\("),pattern = "\\)"),pattern="'"),",")


snp_count = unlist(lapply(alts, function(x){sum(x!="")}))
alts_without_spaces = (lapply(alts, function(x){x[x!=""]}))
vcf_annotation_16$n_alleles = snp_count

ab  =  (lapply(alts_without_spaces, function(x){paste(x,collapse = "_",sep="")}))
vcf_annotation_16$alts = (lapply(alts_without_spaces, function(x){paste(x,collapse = "_",sep="")}))
#vcf_annotation_16$snpid
vcf_in$strand = as.character(strand(vcf_in$X))
vcf_in$corrected_alt = ifelse(vcf_in$strand == "+",vcf_in$varAllele,Biostrings::reverseComplement(vcf_in$varAllele))

vcf_annotation_16$snpid = paste(vcf_annotation_16$chrom,vcf_annotation_16$pos,vcf_annotation_16$alts,sep="_")
vcf_in$snpid = paste(vcf_in$chrom, as.character(vcf_in$start), as.character(vcf_in$corrected_alt),sep="_")
min_prov = unlist(lapply(vcf_in$proveanEffects, function(x){min(x,na.rm=T)}))
min_prov[is.infinite(min_prov)] = NA
vcf_in$min_prov = min_prov
idx_vcf_in = match(vcf_annotation_16$snpid,vcf_in$snpid)

vcf_annotation_16$idx_vcf_in = idx_vcf_in
vcf_annotation_16$min_prov = as.numeric(vcf_in$min_prov[idx_vcf_in])
vcf_annotation_16$joseph_af = as.numeric(vcf_in$maf[idx_vcf_in])
vcf_annotation_16$essential = as.numeric(vcf_in$essential[idx_vcf_in])
vcf_annotation_16 = vcf_annotation_16 %>% mutate(seqnames=chrom,start=pos,end=pos) %>% as_granges()
#### PROVEAN from the other VCF ????? ######
#### Allele-frequecy in the new VCF ????? #### It's probably fine to jjust use the 1002 coding effects file right ###

### Load_all_annotations #### 
### 3004 #### 

skip_binding_sites = T
if(!skip_binding_sites){
yjm981 = list.files("../data/identify_binding_sites/out/YJM981x.fasta/ALIGNED_ENOLOGO_FORMAT_PFMS//",recursive = T,pattern="fimo.tsv",full.names=T)
cbs = list.files("../data/identify_binding_sites/out/CBS2888a.fasta////ALIGNED_ENOLOGO_FORMAT_PFMS//",recursive = T,pattern="fimo.tsv",full.names=T)
by = list.files("../data/identify_binding_sites/out/sacCer3.fasta/ALIGNED_ENOLOGO_FORMAT_PFMS/",recursive  = T, pattern="fimo.tsv",full.names=T)
yps163 = list.files("../data/identify_binding_sites/out/YPS163a.fasta/ALIGNED_ENOLOGO_FORMAT_PFMS/",recursive = T, pattern="fimo.tsv", full.names=T)
yjm454 = list.files("../data/identify_binding_sites/out/YJM454a.fasta/ALIGNED_ENOLOGO_FORMAT_PFMS/",recursive = T, pattern="fimo.tsv", full.names=T)

rm = list.files("../data/identify_binding_sites/out/RMx.fasta//ALIGNED_ENOLOGO_FORMAT_PFMS/",recursive = T, pattern="fimo.tsv", full.names=T)

#length(yjm981)

#cbind(yjm981,cbs,by,yps163, yjm454)

gene_id = unlist(lapply(str_split(yjm981,"/"), function(x){x[9]}))
binding_sites_list = list('3004'=cbind(gene=gene_id,yjm981=yjm981,cbs=cbs),'B'=cbind(gene=gene_id,yps163=yps163,yjm454=yjm454),
                          'A'=cbind(gene=gene_id,by=by,rm=rm))

chromosome_map =  list(chromosome_map=paste("chr",as.roman(1:16),sep=""))

chrom_num_to_roman = data.frame(chrom_num=1:16,chromosome_roman=unlist(chromosome_map))
chrom_num_to_roman = rbind(chrom_num_to_roman,data.frame(chrom_num=17,chromosome_roman="chrM"))



join_overlap_full <- function(x,y){
  options(stringsAsFactors = F)
  hit <- findOverlaps(x,y)
  colnames(mcols(x)) <- paste0(colnames(mcols(x)),".x")
  colnames(mcols(y)) <- paste0(colnames(mcols(y)),".y")
  cns.x <- mcols(x) %>% colnames()
  cns.y <- mcols(y) %>% colnames()
  
  z <- x[queryHits(hit),]
  hitz <- cbind(mcols(y)[subjectHits(hit),],
                mcols(x)[queryHits(hit),])
  colnames(hitz) <- c(cns.y,cns.x)
  mcols(z) <- hitz
  uniqex <- x[-queryHits(hit),]
  uniqey <- y[-subjectHits(hit),]
  mcols(uniqex)[, (length(cns.x)+1):(length(cns.x)+length(cns.y))] <- NA
  colnames(mcols(uniqex)) <- c(cns.x,cns.y)
  mcols(uniqex) <- mcols(uniqex)[,c((length(cns.x)+1):(length(cns.x)+length(cns.y)),1:(length(cns.x)))]
  mcols(uniqey)[, (length(cns.y)+1):(length(cns.x)+length(cns.y))] <- NA
  colnames(mcols(uniqey)) <- c(cns.y,cns.x)
  res <- c(z, uniqex,uniqey)
  mcols(res) <- mcols(res) %>% data.frame() %>% mutate_all(as.character())
  return(res)
}

binding_sites_granges = list()

for (cross in names(binding_sites_list)){
  print(cross)
  cross_df = binding_sites_list[[cross]] %>% as_data_frame()
  for(k in 1:nrow(cross_df)){
    tmp_gene = cross_df$gene[k]
    ds1 = try(
      read.delim(as.character(cross_df[k,2]),sep="\t",comment.char = "#")
    )
    ds2 = try(read.delim(as.character(cross_df[k,3]),sep="\t",comment.char = "#"))
    #ds1
    if(is.null(nrow(ds1))){
      ds1_granges = GRanges()
    }else{
      ds1_filt = ds1 %>% filter(!is.na(stop))
      if (cross == "A"){
        ds1_chrom = ds1_filt$sequence_name
      }else{
        ds1_chrom = chrom_num_to_roman$chromosome_roman[match(ds1_filt$sequence_name,chrom_num_to_roman$chrom_num)]
      }
      #ds1_chrom = chrom_num_to_roman$chromosome_roman[match(ds1_filt$sequence_name, chrom_num_to_roman$chrom_num)]
      ds1_granges = ds1 %>% filter(!is.na(stop)) %>% mutate(seqnames=ds1_chrom,start=start,end=stop)%>% as_granges()
    }
    if(is.null(nrow(ds2))){
      ds2_granges = GRanges()
    }else{
      ds2_filt = ds2 %>% filter(!is.na(stop))
      ## Hack for BY ## 
      ## If
      
      ds2_chrom = chrom_num_to_roman$chromosome_roman[match(ds2_filt$sequence_name,chrom_num_to_roman$chrom_num)]
      ds2_granges = ds2 %>% filter(!is.na(stop)) %>% mutate(seqnames=ds2_chrom, start=start, end=stop) %>% as_granges()
    }
    if(length(ds1_granges) !=0  && length(ds2_granges) !=0){
      outer_m = join_overlap_full(ds1_granges,ds2_granges)
    }else {outer_m = GRanges()}
    binding_sites_granges[[cross]][[tmp_gene]] = list(p1=ds1_granges, p2=ds2_granges,outer_m=outer_m)
  }
}

binding_sites_list_merged =  list()
i = 1

for (cross in names(binding_sites_list)){
  gg = lapply(binding_sites_granges[[cross]], function(x){x$outer_m})
  gg = bind_ranges(gg,.id="NAME")
  parents = crosses.to.parents[[cross]]
  gg$p1 = parents[1]
  gg$p2 = parents[2]
  binding_sites_list_merged[[i]] = gg
  i = i + 1
}
names(binding_sites_list_merged) = names(binding_sites_list)

binding_sites_list_merged[[cross]] %>% filter(q.value.y < 0.05 | q.value.x < 0.05) %>% 
  mutate(score_diff = score.y-score.x ) %>% filter(abs(score_diff) > 2  | (is.na(score.y) | is.na(score.x)))
#gg = lapply(binding_sites_granges[[cross]], function(x){x$outer_m})
}
#bind_ranges(gg,.id="NAME")

#join_overlap_full(binding_sites_granges$A$MAL63_136.meme$p1,binding_sites_granges$A$MAL63_136.meme$p2)
#cbind(yjm=yjm,yps=yps)
#a = read_tsv(yjm[230],comment = "#") %>% mutate(seqnames=sequence_name,start=start,end=stop) %>% as_granges()
#b = read_tsv(cbs[230],comment = "#") %>% mutate(seqnames=sequence_name, start=start,end=stop) %>% as_granges()
#a$score == b$score

#bbb

#g = bbb %>% as_data_frame() %>%  filter(score.y != score.x) %>% mutate(delta_score=abs(score.y - score.x))
#bbb = a %>% join_overlap_left(b) %>% as_data_frame()
#aaa = b %>% join_overlap_left(a) %>% as_data_frame()
#View(rbind(bbb,aaa %>% filter(is.na(matched_sequence.y))) %>% arrange(seqnames, start))
