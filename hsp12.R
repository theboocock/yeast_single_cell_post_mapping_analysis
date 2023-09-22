hsp_locus = data.frame(seqnames="chrVI", start=106166, end=107927 ) %>% as_granges()
cross = "A"
parents =crosses.to.parents[[cross]]
causal_hit= causal_hits%>% filter(shared.parent  %in% parents)  %>% mutate(seqnames=chr) %>% as_granges()
#causal_hit = causal_hit_3004 %>% mutate(seqnames=chr) %>% as_granges()
qtl = qtl_bloom %>% filter(cross == !!cross)
vcf_names = str_remove_all(parents,pattern="a|x")
vcf_annotation = vcf_annotation_16 %>% filter((grepl(vcf_names[1],ref_strains) & grepl(vcf_names[2],alt_strains)) |  (grepl(vcf_names[2],ref_strains) & grepl(vcf_names[1],alt_strains)))
## Missense of STOP
size = unlist(lapply(str_split(vcf_annotation$mutation,"\\|"), function(x){x[3]}))
vcf_annotation$impact  = size 


hsp_locus %>% join_overlap_left(qtl_bloom)
#hsp_locus %>% join_overlap_inner(causal_hits %>% mutate(seqnames=chr) %>% as_granges())
cross_data$A$cis$ASE$ase_cis %>% filter(gene == "YFL014W")
cross_data$A$cis$ASE$ase_noise %>% filter(gene == "YFL014W")
