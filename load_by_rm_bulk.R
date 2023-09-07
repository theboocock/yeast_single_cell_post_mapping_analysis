theta_1000 = readRDS("../rproj/out/combined/Ap/bGLMM_bulk_theta.RDS")
theta_1000_df = data.frame(transcript=names(theta_1000),theta=theta_1000)

gene_expr_1000 = readRDS("../rproj/out/combined/Ap/mean_TPM_bulk_BYxRM1000.RDS")
gene_expr_1000_df =data.frame(transcript=names(gene_expr_1000),expr=gene_expr_1000)
cis_test_1000 = readRDS("../data/bulkEQTL_cis_only_test.RDS")
cis_test_1000$p_adj = p.adjust(cis_test_1000$p.value,method="fdr")
#(("../rproj/out/combined/Ap/"))
bulk_local_eQTL = cis_test_1000%>% inner_join(gene_expr_1000_df,by="transcript")
augmented_peaks = readRDS("../rproj/out/peaks_per_gene_augmented.RDS")
aug_peaks = bind_rows(augmented_peaks,.id="gene")

get_chrom_pos_end = function(x){
  x2 = str_split(x,":")
  chrom = unlist(lapply(x2, function(x){x[1]}))
  pos_str = unlist(lapply(x2, function(x){x[2]}))
  pos = unlist(sapply(str_split(pos_str,"_"), function(x){x[1]}))
  return(data.frame(chrom,pos))  
}


cross_data$A_bulk$cis = bulk_local_eQTL

segdata= readRDS("../rproj/out/combined/A/segData.RDS")
pos_chr_peak = get_chrom_pos_end(aug_peaks$pmarker)
gff_in_df = gff_in %>% as_data_frame()
gff_in_df = gff_in_df %>% mutate(mid=round((start+end)/2))
gff_in_df  = gff_in_df[gff_in_df$seqnames != "chrmt",]
gff_in_df$seqnames = factor(gff_in_df$seqnames)
aug_peaks$mid = gff_in_df$mid[match(aug_peaks$gene,gff_in_df$Name)]
aug_peaks$tchr = gff_in_df$seqnames[match(aug_peaks$gene,gff_in_df$Name)]
aug_peaks$chrom = pos_chr_peak$chrom
aug_peaks$pos = as.numeric(pos_chr_peak$pos)
aug_peaks$tchr = factor(aug_peaks$tchr, rev(levels(aug_peaks$tchr)))
aug_peaks$chrom = factor(aug_peaks$chrom,levels=levels(seqnames(gff_in))[-17])

chrom = unlist(lapply(str_split(colnames(segdata$Gsub),"_"), function(x){x[1]}))
pos = unlist(lapply(str_split(colnames(segdata$Gsub),"_"), function(x){x[2]}))
markerGR = data.frame(seqnames=chrom,start=as.numeric(pos),width=1) %>% as_granges()
cmin=sapply(sapply(split(markerGR, seqnames(markerGR)), start), min)
cmin[!is.na(cmin)]=0
cmax=sapply(sapply(split(markerGR, seqnames(markerGR)), start), max)
cbin=data.frame(chr=names(cmin),
                start=cmin,
                end=cmax, 
                strand="*")
cbin=makeGRangesFromDataFrame(cbin)
cbin50k=getcbin50k_genome()
cPf = aug_peaks %>% filter(tchr != chrom)

bin.table = makeBinTable(cPf, cbin50k)
sig.hp=qpois(1-(.05/length(bin.table$pos)),ceiling(mean(bin.table$count)))+1
bin.table %>% filter(count > sig.hp)
cPf=getBinAssignment(cPf, cbin50k,markerGR )
sig.hp.names=table(cPf$bin)>sig.hp
cPf$in.hotspot=sig.hp.names[cPf$bin]
#cbin_table = makeBinTable(cPf,cbin50k)
cPf$chrom_short = str_replace_all(cPf$chr,"chr","")
cPf$tchrom_short = str_replace_all(cPf$tchr,"chr","")
cPf$transcript = cPf$gene
cPf$cell_cycle_stage = "combined"
cPf$tpos = cPf$gene.gcoord
order_levels = c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")
rev_order_levels = rev(order_levels)
cPf$chrom_short_f = factor(cPf$chrom_short,levels=order_levels)
cPf$tchrom_short_f = factor(cPf$tchrom_short,levels=rev_order_levels)

bulk_hotspots = cPf
cross = "A"

### Something custom for this cross ###
#bulk_hotspots$transcript = bulk_hotspots$gene
#bulk_hotspots$cell_cycle_stage= "combined"

cil = unlist(as.numeric(lapply(str_split(bulk_hotspots$CI.l,"_"), function(x){str_split(x[1],":")[[1]][2]})))
cir =  unlist(as.numeric(lapply(str_split(bulk_hotspots$CI.r,"_"), function(x){str_split(x[1],":")[[1]][2]})))
bulk_hotspots$CI.l = cil
bulk_hotspots$CI.r = cir
bulk_hotspots$peak.marker = bulk_hotspots$pmarker
bulk_hotspots$peak.marker = str_replace_all(bulk_hotspots$peak.marker,":","_")

#background = (names(gene_expr_1000))

#lod_df_beta, cell_cycle_beta_df2,segdata,out_name

hotspot_peaks = bulk_hotspots
combined_peaks = cPf
lod_df_beta = cross_data$A$trans$cell_cycle_lods
cell_cycle_beta_df2 = cross_data$A$trans$cell_cycle_beta_all
segdata=segdata
background = names(gene_expr_1000)
out_name = "A_old"


hotspots_bulk = hotspot_enrichment_and_function(cross="A",hotspot_peaks = bulk_hotspots,combined_peaks = cPf,lod_df_beta = cross_data$A$trans$cell_cycle_lods,
                                                cell_cycle_beta_df2=cross_data$A$trans$cell_cycle_beta_all,segdata = segdata, background = names(gene_expr_1000),fdr_fx = NULL,out_name = "A_old")




hotspot_bins_final = hotspots_bulk$hotspot_list %>% group_by(bin,hotspot_pos, chrom) %>% summarise(n=n())
#hotspot_bins_final = enrich_list$hotspot_list %>% group_by(bin,hotspot_pos,chrom) %>% summarise(n=n())

hotspot_bins_final$chrom = unlist(lapply(str_split(hotspot_bins_final$bin,":"), function(x){x[1]}))
pos_str = unlist(lapply(str_split(hotspot_bins_final$bin,":"), function(x){x[2]}))
start = as.numeric(unlist(lapply(str_split(pos_str,"-"), function(x){x[1]})))
end =as.numeric(unlist(lapply(str_split(pos_str,"-"), function(x){x[2]})))
hotspot_bins_final$start = start 
hotspot_bins_final$end = end
hotspot_bins_final_gr = makeGRangesFromDataFrame(hotspot_bins_final)
hotspot_bins_final_gr$n = hotspot_bins_final$n
hotspot_bins_final_gr$bin = hotspot_bins_final$bin
hotspot_bins_final_gr$hotspot_pos = hotspot_bins_final$hotspot_pos


#bin.table

cross_data$A_bulk$cis = bulk_local_eQTL
cross_data$A_bulk$trans = list()
cross_data$A_bulk$trans[["combined_peaks"]] = bulk_hotspots
cross_data$A_bulk$trans[["hotspot_threshold"]] = sig.hp
cross_data$A_bulk$trans[["hotspot_peaks"]]= bulk_hotspots %>% filter(in.hotspot)
cross_data$A_bulk$trans[["segdata"]] = segdata
cross_data$A_bulk$trans[["hotspot_list"]] = hotspots_bulk$hotspot_list
cross_data$A_bulk$trans[["hotspot_enrichments_and_overlaps"]] = hotspots_bulk$annotation_list
cross_data$A_bulk$trans[["hotspot_table"]] = bin.table
cross_data$A_bulk$trans[["hotspot_table"]]$chrom_f = convert_chrom_to_simple_factor(cross_data$A_bulk$trans[["hotspot_table"]]$chr)


cross_data$A_bulk$trans[["hotspot_bins_final"]] = hotspot_bins_final_gr



#cross="A"
#hotspot_peaks = bulk_hotspots
#lod_df_beta = cross_data$A$trans$cell_cycle_lods
#cell_cycle_beta_df2=cross_data$A$trans$cell_cycle_beta_all
#segdata = segdata
#out_name = "A_old"
(table(hotspots_bulk$hotspot_list$bin))
length(table(cross_data$A$trans$hotspot_peaks$bin[cross_data$A$trans$hotspot_peaks$in.hotspot]))

(table(cross_data$A$trans$hotspot_list$bin))

cross_data$A$trans$hotspot_table %>% ggplot(aes(y=count,x=as.factor(pos))) + geom_col()  + geom_hline(yintercept = cross_data$A$trans$hotspot_threshold)+  facet_wrap(~chr,scales="free")

bin.table %>% ggplot(aes(y=count,x=as.factor(pos))) + geom_col()  + geom_hline(yintercept = sig.hp)+  facet_wrap(~chr,scales="free")
