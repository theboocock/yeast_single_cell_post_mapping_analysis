source("plate_reader_utils.R")
conditions = "/media/theboocock/Data/Dropbox/PHDTHESIS/projects/gal_final_github_october/data/plate_reader/mutants_final_freezer_gal/gal_mutants_glu_feb2019_conditions.csv"
gpa_validation_folder = "/media/theboocock/Data/Dropbox/Postdoc/projects/gpa_validation/"
glue("{gpa_validation_folder}/")


add_labels_laura_strains2 = function(ypd_growth){
  ypd_growth$df$ind2=as.character(ypd_growth$df$ind)
  ypd_growth$df$ind2[ypd_growth$df$ind == "416_GLU"] = "WT"
  ypd_growth$df$ind2[ypd_growth$df$ind == "417_GLU"] = "WT"
  ypd_growth$df$ind2[ypd_growth$df$ind == "420_GLU"] = "82R"
  ypd_growth$df$ind2[ypd_growth$df$ind == "421_GLU"] = "82R"
  ypd_growth$df$ind2[ypd_growth$df$ind == "422_GLU"] = "82R"
  ypd_growth$df$ind2[ypd_growth$df$ind == "423_GLU"] = "82R"
  ypd_growth$df$ind2[ypd_growth$df$ind == "424_GLU"] = "469I"
  ypd_growth$df$ind2[ypd_growth$df$ind == "425_GLU"] = "469I"
  return(ypd_growth)
}

add_labels_laura_strains = function(ypd_growth){
  ypd_growth$df$ind2=as.character(ypd_growth$df$ind)
  print(ypd_growth$df$ind2 == "1_GLU")
  ypd_growth$df$ind2[ypd_growth$df$ind == "1_GLU"] = "469I"
  ypd_growth$df$ind2[ypd_growth$df$ind == "2_GLU"] = "WT"
  ypd_growth$df$ind2[ypd_growth$df$ind == "3_GLU"] = "469I 82R"
  ypd_growth$df$ind2[ypd_growth$df$ind == "4_GLU"] = "82R"
  ypd_growth$df$ind2[ypd_growth$df$ind == "5_GLU"] = "82R"
  ypd_growth$df$ind2[ypd_growth$df$ind == "6_GLU"] = "469I 82R"
  print(ypd_growth$df$ind2)
  
  return(ypd_growth)
}



samples1 = glue("{gpa_validation_folder}/data/plate_leslie_july23/samples.csv")
#mutant_pheno=  "/media/theboocock/Data/Dropbox/PHDTHESIS/projects/gal_final_github_october/ext/data/gal_mutant_pheno_nov.csv" 

plate1 = glue("{gpa_validation_folder}data/plate_leslie_july23/072423_Plate1.csv")
plate2 = glue("{gpa_validation_folder}data/plate_leslie_july23/072423_Plate2.csv")
plate3 = glue("{gpa_validation_folder}data/plate_leslie_july28/plate1.csv")
plate4 = glue("{gpa_validation_folder}data/plate_leslie_july28/plate2.csv")

ypd_growth1 = add_labels_laura_strains2(process_phenotype_plates(plate1, conditions, samples1))
ypd_growth2 = add_labels_laura_strains2(process_phenotype_plates(plate2, conditions, samples1))

ypd_growth3 = add_labels_laura_strains2(process_phenotype_plates(plate3, conditions, samples1))
ypd_growth4 = add_labels_laura_strains2(process_phenotype_plates(plate4, conditions, samples1))



samples2= glue("{gpa_validation_folder}data/lwm_samples.csv")

plate5 = glue("{gpa_validation_folder}data/ypd_1_13.csv")
plate6 = glue("{gpa_validation_folder}data/ypd_1_18.csv")
plate7 = glue("{gpa_validation_folder}data/ypd_1_23_lidonly.csv")


ypd_growth5 = add_labels_laura_strains(process_phenotype_plates(plate5, conditions, samples2))
ypd_growth6 = add_labels_laura_strains(process_phenotype_plates(plate6, conditions, samples2))
ypd_growth7 = add_labels_laura_strains(process_phenotype_plates(plate7, conditions, samples2))

plate5 = ypd_growth5$df %>% dplyr::select(ind2,values, doubling) %>%  dplyr::rename(ID=ind2)
plate5$plate = "Laura 1_13_23"
plate6 = ypd_growth6$df %>% dplyr::select(ind2,values, doubling)%>%  dplyr::rename(ID=ind2)
plate6$plate = "Laura 1_18_23"
plate7 = ypd_growth7$df %>% dplyr::select(ind2,values, doubling)%>%  dplyr::rename(ID=ind2)
plate7$plate= "Laura 1_23_23"
plate1 = ypd_growth1$df  %>% filter(values < 0.085 & values > 0.07)%>% dplyr::select(ind2,values,doubling) %>% dplyr::rename(ID=ind2)
plate1$plate = "Leslie 1 7_26_23"
plate2 =  ypd_growth2$df %>% filter(values < 0.085 & values > 0.07)%>% dplyr::select(ind2,values,doubling) %>% dplyr::rename(ID=ind2)
plate2$plate = "Leslie 2 7_26_23"
plate3 = ypd_growth3$df%>% filter(values < 0.085 & values > 0.07) %>%dplyr::select(ind2,values,doubling) %>% dplyr::rename(ID=ind2)
plate3$plate = "Leslie 1 7_28_23"
plate4 = ypd_growth4$df%>% filter(values < 0.085 & values > 0.07)%>% dplyr::select(ind2,values,doubling) %>% dplyr::rename(ID=ind2)
plate4$plate = "Leslie 2 7_28_23"
#plate4

p_m = rbind(plate1,plate2,plate3,plate4,plate5,plate6,plate7) #%>% ggplot(aes(y=))
p_m = p_m %>% filter(ID != "469I 82R")
p_m %>% filter(values < 0.085 & values > 0.07)
m_resid = residuals(lm((p_m$doubling) ~ p_m$plate))
p_m$resid = m_resid + 0.4519579

mycomparsions = c("WT","W82R")
p_m$strain_new = factor(p_m$ID, levels=c("WT","82R","469I"))

mycomparsions = list(c("WT","82R"),c("82R","469I"),c("WT","469I"))

pa = p_m %>% mutate(plate_l = grepl("Les",plate)) %>% 
  ggplot(aes(y=resid,x=strain_new)) + geom_boxplot() + geom_jitter(width=0.2) + theme_bw() + ylab("Doublings per hour") + xlab("Strain")  + theme(text=element_text(size=18)) +
  stat_compare_means(comparisons = mycomparsions,method="t.test",size=8) 
m1 = (lm((p_m$doubling) ~ p_m$ID + p_m$plate))

library(emmeans)
em = emmeans::emmeans(m1,specs=~ID)
pairs(em)
#p#ypd_growth1$


mating1 = read.csv("data/mating1.csv")
mating2=read.csv("data/mating2.csv")


mating1$exp = 1
mating2$exp = 2
mating_df = rbind(mating1,mating2)
s1 = unlist(lapply(str_split(mating_df$Strain," "), function(x){x[1]}))
s2 = unlist(lapply(str_split(mating_df$Strain," "), function(x){x[2]}))
blue = grepl("B",s1)
red = grepl("R",s1)

s1 = str_replace_all(s1,"R|B","")
s2 = str_replace_all(s2,"R|B","")


mating_df= mating_df%>% mutate(blue_strain = ifelse(blue,s1,s2))
mating_df=  mating_df %>% mutate(red_strain = ifelse(red,s1,s2))

mating_df$strain2 = factor(mating_df$Strain,levels=c("416B 420R","416R 420B","416B 424R","416R 424B","420B 424R","420R 424B"))

mating_df$expr = c(rep(1,12),rep(2,12),rep(3,12),rep(4,12),rep(5,12))
aa  = data.frame()
mating_df$XX = residuals(lm(mating_df$X ~ mating_df$color))
for(j in unique(mating_df$expr)){
  print(j)
  b = mating_df %>% filter(expr == j & color=="red")
  norm = mean(b$XX[which(b$red_strain == "416")])
  b$norm = b$XX - norm + 100
  b$strain = b$red_strain
  b$color2 = "red"
  aa = rbind(aa,b)
  b = mating_df %>% filter(expr == j & color=="blue")
  norm = mean(b$XX[which(b$blue_strain == "416")])
  b$norm = b$XX - norm + 100
  b$strain = b$blue_strain
  b$color2 = "blue"
  aa = rbind(aa,b)
}

aa = aa %>% mutate(strain_new = case_when(
  strain =="416" ~ "WT",
  strain == "420" ~ "82R",
  strain == "424" ~ "469I"))
aa$strain_new = factor(aa$strain_new, levels=c("WT","82R","469I"))
pb = aa%>% ggplot(aes(y=norm,x=strain_new)) + geom_boxplot(outlier.shape  = NA) +  geom_jitter(width=.2,size=3) +
  ylab("Mating efficiency") + xlab("Strain") + scale_color_manual(values=c("blue","red")) + theme_bw()  +
  theme(text=element_text(size=18))  +stat_compare_means(comparisons = mycomparsions,method="t.test",size=8)#+ ylim(c(99,112))


aa_red = aa %>% filter(color=="red")
aa_blue = aa %>% filter(color=="blue")
mm = lm(aa$norm ~ aa$strain_new)
aaaa = emmeans(mm,spec=~strain_new)
pairs(aaaa)

plot_side_left = pa / pb
plot_side_left
ggsave("figures/figure_6_gpa.svg",width=8,height=12)
library(SNPRelate)
genofile = snpgdsOpen("1012.gds")
#snpgdsClose(genofile)
dissim = snpgdsDiss(genofile,autosome.only = F)
rownames(dissim$diss)= dissim$sample.id
colnames(dissim$diss) = dissim$sample.id
joseph_annotation = xlsx::read.xlsx("/media/theboocock/Data/Dropbox/Postdoc/projects/all_yeast_alignment/align_avaliable_yeast_genomes_2022/full_annotations.xlsx",sheetName="peter2018")

joseph_annotation$Clades_trim = (str_trim(joseph_annotation$Clades))
joseph_annotation$mosaic = grepl("Mosaic",joseph_annotation$Clades)
m1 = joseph_annotation[joseph_annotation$mosaic,]

###
idx_mosaic = which(rownames(dissim$diss) %in% m1$standardized_name2)
tree = ape::bionj(dissim$diss)
mid = phangorn::midpoint(tree)
pos = read.gdsn(index.gdsn(genofile,"snp.position"))
chrom =  read.gdsn(index.gdsn(genofile,"snp.chromosome"))
allele =  read.gdsn(index.gdsn(genofile,"snp.allele"))

#Gpa1
idx = which(chrom == "VIII" & pos == "114674")
#idx2 = which(chrom == "XIII" & pos == "25025")


gt =  read.gdsn(index.gdsn(genofile,"genotype"))
w82r_gt = (2-(gt[,idx]))
gt[gt==3] = NA
gt3 = 2-gt
ac = colSums(gt3,na.rm=T)
tc = colSums(!is.na(gt3))
af = ac / (tc*2)
gt3[,af > 0.5] = 2- gt3[,af  > .5]
maf = ifelse(af > 0.5,1-af,af)
mosaic_prop = sum(w82r_gt[idx_mosaic])/(length(idx_mosaic) * 2 )

het_sites = rowSums(gt == 1,na.rm=T)
het_sites_divisor = rowSums(!is.na(gt))*2
het_sites_prop = het_sites/het_sites_divisor

het_df = data.frame(het_sites_prop=het_sites_prop,gpa=w82r_gt)
het_df$sample = dissim$sample.id
het_df = het_df %>% inner_join(joseph_annotation,by=c("sample"="standardized_name2"))


ccc = het_df %>% filter(gpa != 1) %>% filter(Ploidy != 1)
fisher.test(table(ccc$gpa,ccc$Zygosity))

anova(lm(ccc$het_sites_prop ~ ccc$gpa))

wilcox.test(ccc$het_sites_prop ~ ccc$gpa)

props = c()
for(idx2 in sample(idx_a/ll,size=1e4,replace = T)){
  tmp_gt2 = gt3[,idx2]
  #i  
  len =  sum(!is.na(tmp_gt2[idx_mosaic]))
  mosaic_prop2 = sum(tmp_gt2[idx_mosaic],na.rm=T)/(len * 2 )
  #mosaic_prop2= ifelse(mosaic_prop2 > 0.5,1-mosaic_prop2,mosaic_prop2)
  props = c(props, mosaic_prop2)
}

clade_summary = het_df %>% group_by(Clades) %>% summarise(maf=sum(gpa)/(n()*2),n=n())
library(ggtree)
library(treeio)
aa_2 = as.treedata(mid)
gfive = groupOTU(aa_2@phylo, list(het_df$sample[het_df$gpa == 2],het_df$sample[het_df$gpa == 1],het_df$sample[het_df$gpa == 0]))
p3 = ggtree(gfive, aes(color=group), layout="fan",open.angle = 180 )  + scale_color_manual(values=c("#377eb8","#e41a1c","grey85"))
ggsave("figures/tree.svg")

clade_list = list()
clade_list_new_names = list()
i = 1
new_clades = unlist(lapply(str_split(unique(het_df$Clades),"\\."),function(x){x[2]}))
new_clades[is.na(new_clades)]= "Missing"
x = 1
for(clade in unique(het_df$Clades)){
  print(clade)
  clade_list[[i]] = het_df$sample[het_df$Clades == clade]
  clade_list_new_names[[i]] = het_df$sample[het_df$Clades == clade]
  i = i + 1
}
names(clade_list) = unique(het_df$Clades)
names(clade_list_new_names) = new_clades

gfour = groupOTU(aa_2@phylo, clade_list)

#clade_summary %>% 
p4 = clade_summary %>% ggplot(aes(x=maf,y=reorder(Clades,maf))) + geom_bar(stat="identity")  + theme_classic() + theme(text=element_text(size=18)) + xlab("82R Allele-frequency") + ylab("Clades")

ggtree(gfour,aes(color=group),layout="fan")


((pa / pb ) | (p3/p4)) + plot_annotation(tag_levels = 'A')

ggsave("figure5.png")
#(pa / pb ) | p3


