library(SNPRelate)
library(tidyverse)

#snpgdsVCF2GDS("/media/theboocock/scratch/vcf_tmp/maf_0.01_joseph.vcf.gz","rare.gds",method="biallelic.only")
genofile2 = snpgdsOpen("rare.gds") # MAF > 0.01 & MAF < 0.04??
pos2 = read.gdsn(index.gdsn(genofile2,"snp.position"))
chrom2 =  read.gdsn(index.gdsn(genofile2,"snp.chromosome"))
allele2 =  read.gdsn(index.gdsn(genofile2,"snp.allele"))

idx_rare = which(pos2 == 113512 & chrom2 == "VIII")
gt2 =  read.gdsn(index.gdsn(genofile2,"genotype"))
gt_by = gt2[,idx_rare]
#i#dx_mosaic2 = which(rownames(dissim$diss) %in% m1$standardized_name2)
snpgdsVCF2GDS("/media/theboocock/scratch/vcf_tmp/maf_0.05_joseph.vcf.gz","1012.gds",method="biallelic.only")
snpgdsSummary("1012.gds")
genofile = snpgdsOpen("1012.gds")
#snpgdsClose(genofile)
dissim = snpgdsDiss(genofile,autosome.only = F)
rownames(dissim$diss)= dissim$sample.id
colnames(dissim$diss) = dissim$sample.id

### 


joseph_annotation = xlsx::read.xlsx("/media/theboocock/Data/Dropbox/Postdoc/projects/all_yeast_alignment/align_avaliable_yeast_genomes_2022/full_annotations.xlsx",sheetName="peter2018")

joseph_annotation$mosaic = grepl("Mosaic",joseph_annotation$Clades)
#joseph_annotation

joseph_annotation$Geographical.origins
colnames(joseph_annotation)
m1 = joseph_annotation[grep("Mosaic",joseph_annotation$Clades),]
m1$standardized_name2
by_var  = data.frame(sample=rownames(dissim$diss),gt=gt_by) %>%
  inner_join(joseph_annotation,by=c("sample"="standardized_name2"))

by_var %>% group_by(Clades) %>% summarise(gt_freq=sum(gt)/(2*n()),n=n())  %>% ggplot(aes(y=Clades,x=gt_freq)) +
  geom_bar(stat="identity")

#gt_tmp

## Rare variant first #
# = joseph_annotation$Clades[grep("Mosaic",joseph_annotation$Clades)]

idx_mosaic = which(rownames(dissim$diss) %in% m1$standardized_name2)
idx_mosaic_num = ifelse(rownames(dissim$diss) %in% m1$standardized_name2,1,0)

non_hap = by_var %>% filter(Ploidy != 1)
fisher.test(table(non_hap$Zygosity[non_hap$gt != 1],non_hap$gt[non_hap$gt!=1]))

fisher.test(table(by_var$Zygosity[by_var$gt != 1],by_var$gt[by_var$gt != 1]))

View(by_var)

table(by_var$Zygosity[by_var$gt != 1],by_var$gt[by_var$gt != 1])
#glm(by_var$Zygosity ~ by_var$gt,family="binomial")

#idx_mosaic
by_var %>% ggplot(aes(y=Proportion.of.clean.heterozygous.SNPs..whole.dataset.,x=factor(gt))) + geom_boxplot()


tree = ape::bionj(dissim$diss)
mid = phangorn::midpoint(tree)
pos = read.gdsn(index.gdsn(genofile,"snp.position"))
chrom =  read.gdsn(index.gdsn(genofile,"snp.chromosome"))
allele =  read.gdsn(index.gdsn(genofile,"snp.allele"))

idx = which(chrom == "VIII" & pos == "114674")

gt =  read.gdsn(index.gdsn(genofile,"genotype"))
tmp_gt = (2-(gt[,idx]))
gt[gt==3] = NA
gt3 = 2-gt





#which(chrom == "VIII" &	pos == "113512")

ab = colSums(gt3,na.rm = T)

idxs = which(ab > 385 & ab < 415 + 20 )
idxs = which(ab > 385 & ab < 415 + 20 )
ab = -log10((pnorm(-abs(zs)))*2)
ac = colSums(gt3,na.rm=T)
tc = colSums(!is.na(gt3))
af = ac / (tc*2)
maf = ifelse(af > 0.5,1-af,af)


gt3[,af > 0.5] = 2- gt3[,af  > .5]

sum(colSums(gt3,na.rm=T) > nrow(gt3))
idx_all = (maf < maf[idx] + .02 & maf > maf[idx] -0.02)
idx_all[idx] =F
idx_all = which(idx_all)
#maf >

mosaic_prop = sum(tmp_gt[idx_mosaic])/(length(idx_mosaic) * 2 )

mosaic_prop_by= sum(gt_by[idx_mosaic])/(length(idx_mosaic) * 2 )

mosaic_maf = 1-mosaic_prop

het_sites = rowSums(gt == 1,na.rm=T)
het_sites_divisor = rowSums(!is.na(gt))*2
het_sites_prop = het_sites/het_sites_divisor
#hist(het_sites/het_sites_divisor)
cor(tmp_gt,het_sites_prop,method="spearman")

summary(lm(scale(het_sites_prop) ~ (tmp_gt)))

data.frame(het_sites_prop=het_sites_prop,gpa=tmp_gt) %>% ggplot(aes(y=het_sites_prop,x=factor(gpa))) + geom_boxplot() + geom_jitter()


aaa = data.frame(sample=colnames(dissim$diss), het_sites_prop=het_sites_prop,gpa=tmp_gt) #%>% filter(gpa == 2)# %>% het_sites_prop # %>% 

bbb = (aaa %>% inner_join(joseph_annotation,by=c("sample"="standardized_name2"))) #%>% filter()

View(bbb %>% group_by(Clades) %>% summarise(n=n()))

View(bbb %>% group_by(Clades) %>% summarise(n=sum(gpa)/(n()*2)))
bbb = bbb %>% mutate(mosaic_mixed=sample %in% m1$standardized_name2)

(bbb %>% group_by(mosaic_mixed) %>% summarise(n=sum(gpa)/(n()*2)))


bbb %>% ggplot(aes(x=factor(gpa),y=Zygosity)) + geom_point()
ccc = bbb %>% filter(gpa != 1) %>% filter(Ploidy != 1)
fisher.test(table(ccc$gpa,ccc$Zygosity))
#fisher.test(table(bbb$gpa,bbb$Zygosity))
#bbb %>% 

props = c()
for(idx2 in sample(idx_all,size=1e4,replace = T)){
  tmp_gt2 = gt3[,idx2]
#i  
  len =  sum(!is.na(tmp_gt2[idx_mosaic]))
  mosaic_prop2 = sum(tmp_gt2[idx_mosaic],na.rm=T)/(len * 2 )
  #mosaic_prop2= ifelse(mosaic_prop2 > 0.5,1-mosaic_prop2,mosaic_prop2)
  props = c(props, mosaic_prop2)
}

#mosaic_maf
sum(mosaic_prop < props)

#library(tidyverse)
data.frame(props=props) %>% ggplot(aes(x=props)) + geom_histogram() + geom_vline(xintercept = mosaic_prop)

#mmp model.matrix(lm(idx_mosaic_num~tmp_gt2))
library(fastglm)
model.matrix(tmp_gt2)
#for(i in 1:)
zs = c()
beta= c()
se=c()
coefs = data.frame()
for(i in 1:ncol(gt)){
  print(i)
  tmp_gt2 = (2-gt[,i])
  idx_na = is.na(tmp_gt2)
  mmp = model.matrix(lm(idx_mosaic_num[!idx_na]~scale(tmp_gt2[!idx_na])))
  
  assoc = (fastglm(y=idx_mosaic_num[!idx_na],x=mmp,family="binomial"))
  
  m3 = summary(lm(het_sites_prop[!idx_na] ~ (tmp_gt2[!idx_na])))
  coef = data.frame(t(m3$coefficients[2,1:2]))
  coefs = rbind(coefs,coef)
  z = assoc$coefficients[2]/assoc$se[2]
  zs = c(zs,z)
  beta = c(beta, assoc$coefficients[2])
  se = c(se,assoc$se[2])
  
}

ab = -log10((pnorm(-abs(zs)))*2)
ac = colSums(gt3,na.rm=T)
tc = colSums(!is.na(gt3))
af = ac / (tc*2)
maf = ifelse(af > 0.5,1-af,af)

data.frame(beta=beta,se=se,zs=zs,chrom=chrom,pos=pos,log10p=ab,af=af) %>% filter(beta > -2.5) %>% ggplot(aes(y=log10p,x=pos))  + geom_point()+ 
  facet_wrap(~chrom,scales="free")

coefs[idx,]
coefs  %>% mutate(chrom=chrom,pos=pos) %>% mutate(z=Estimate/Std..Error) %>% filter(abs(z) > 20)
res  = data.frame(beta=beta,se=se,zs=zs,chrom=chrom,pos=pos,log10p=ab,af=af,maf=maf) 
res %>% filter(maf > .15) %>% filter(abs(zs) > 8) %>% arrange(-log10p)



gt_gpa1 = gt[,34395]

aa = cor(gt_gpa1,gt[,30000:40000])



gt_br = colnames(dissim$diss) %in% joseph_annotation$standardized_name2[grep("Bra",joseph_annotation$Clades)]
aaa = gt[gt_br,]


aaa

plot(aaa[1,])

ac = colSums(aaa,na.rm=T)
tc = colSums(!is.na(aaa))
keep_vars = !(ac/(2*tc) > 0.9 | ac/(2*tc)<0.1)
plot(aaa[3,keep_vars])

rowSums(aaa==1,na.rm=T)

data.frame(chrom=chrom[keep_vars],pos=as.numeric(pos[keep_vars]),gt=aaa[8,keep_vars]) %>% ggplot(aes(y=gt,x=pos)) + 
  facet_wrap(~chrom,scales="free") + geom_point()

library(milorGWAS)
x <- read.bed.matrix(basename = "/media/theboocock/scratch/vcf_tmp/maf_0.05")
x1 =  read.bed.matrix(basename = "/media/theboocock/scratch/vcf_tmp/maf_0.01")

K = GRM(x)

ordered_pheno = joseph_annotation[match(x@ped$id,joseph_annotation$standardized_name2,),]


x@ped$id
ordered_pheno$standardized_name2
aa = as.bed.matrix(x)
assoc_test = milorGWAS::association.test.logistic(x=x1,
                                     Y=ordered_pheno$mosaic, 
                                     X=matrix(1,length(ordered_pheno$mosaic),1),K=K, algorithm='offset')


go = assoc_test[order(assoc_test$p),]%>% filter(chr == "chrVIII" & (pos > 112512 & pos < 115600))# & pos == 114674)


assoc_test[order(assoc_test$p),]%>% filter(p < 1e-18)

#plot(-log10(assoc_test$p))
aa$X$prov_min = aa$prov_min
bb = aa$X %>% as_data_frame() 


al = bb %>% inner_join(go,by=c("start"="pos"))
al2 = al[al$CONSEQUENCE  == "nonsynonymous",]
al2$prov_min
View(al2)
