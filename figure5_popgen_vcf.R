coding = readRDS("/media/theboocock/scratch/vcf_tmp/vcf_codingEffects.RDS")
library(Biostrings)
prov_min_all  = unlist(lapply(coding$proveanEffects, min))
coding$prov_min = prov_min_all
coding2 = coding[!is.na(prov_min_all),]

gl = readRDS("/media/theboocock/scratch/vcf_tmp/vcf_granges.RDS")
prov_min = unlist(lapply(aa$proveanEffects, function(x){min(x)}))
#
#aa

aa = coding[coding$GENEID == "YHR005C",]
bb = coding[coding$GENEID == "YKR084C",]
af = gl[match(names(bb$X),names(gl$X)),]

gl

coding$CONSEQUENCE

aa$prov = prov_min
prov = aa[unlist(nchar(aa$REF) == 1 & nchar(aa$ALT) == 1),]
library(plyranges)i
prov %>% filter(CONSEQUENCE == "nonsynonymous")
bb = prov[prov$CONSEQUENCE == "nonsynonymous",]
plot(y=bb$prov , x=bb$PROTEINLOC)



gl$ chrXIII:25025_A/G


coding$X[which(coding$X$GENEID == "YML123C"),] %>% filter(REFAA == "L") %>% filter(VARAA=="P")
#bb$

sum(duplicated(coding$X))
bb2 = bb[!duplicated(bb$X),]
bb2$prov
bb2

which(names(gl) %in% names(bb2))
af = gl[match(names(bb$X),names(gl$X)),]
plot(bb$prov,af$maf)


af$maf > 0.1

bb[af$maf > .1,]


af = gl$af[match(names(coding2$X), names(gl$X))]
maf = gl$maf[match(names(coding2$X), names(gl$X))]


#gl
#coding
plot(coding2$prov_min,maf)


#af

coding2$prov_min < -5

coding2[coding2$prov_min < -10 & maf > .15,]
