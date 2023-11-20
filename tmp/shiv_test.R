library(ShortRead)


a = readFastq("~/Downloads/T2GAL1_S2_L001_R1_001.fastq")

ab =  sread(a)


ab[1:6,]

mean(table(substr(ab,start=1,stop = 6)))


mean(table(substr(ab,start=1,stop = 6)))


guides = (substr(ab,start=27,stop=27+19))
names(tab)[min(tab) == tab]

