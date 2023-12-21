

strains = unique(c("YJM451","YJM1399","YJM1386","YJM554","YJM1339","YJM1418","YJM1252","YJM1444","YJM1338","YJM451"))

grep_str= paste(strains,sep="",collapse = "|")


het_df[grep(grep_str,het_df$Isolate.name),]
