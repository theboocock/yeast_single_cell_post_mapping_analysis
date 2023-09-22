

ap_match = readRDS("../rproj/out/combined/Ap/segMatchList.RDS")
ap_match_round =  ifelse(ap_match$vg> 0.5,1,0)
ap_match_round2 = ifelse(ap_match$sprevG > 0.5,1,0)

ap_match_combined = t(ap_match_round2[,ap_match$best_match_seg])
#ap_match-
#i#ap_match_combined

accuracy = apply(ap_match_combined == ap_match_round, 1, sum)/ncol(ap_match_combined)
ap_df = readRDS("../rproj/out/combined/Ap/segData.RDS")
ap_df$barcode.features$accuracy = accuracy

ap_df$barcode.features %>% ggplot(aes(y=accuracy * 100,x=(nUMI))) + geom_point() +  scale_x_log10() + xlab(expression('UMI count')) + ylab("Genotyping Accuracy") + theme_bw() + stat_cor(method="spearman",size=12,cor.coef.name = "rho") +
  theme(text=element_text(size=24)) + ylim(c(50,100)) + geom_hline(yintercept = median(ap_df$barcode.features$accuracy*100),color="red")
ggsave("fig_final/s5.png",width=16,height = 12,dpi=300)
