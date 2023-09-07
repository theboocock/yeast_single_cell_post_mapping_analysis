gsub=cross_data$`3004`$trans$segdata$Gsub
idx_one = c(which(colnames(gsub) == "chrXI_63984_C_T_45553") -1, which(colnames(gsub) == "chrXI_63984_C_T_45553") +1)
colnames(gsub)[idx_one]
