seg_match = readRDS("../rproj/out/combined/Ap/segMatch.RDS")
seg_match_df = data.frame(seg=names(table(seg_match$best_match_seg)),table(seg_match$best_match_seg))
seg_match_df %>% ggplot(aes(x=Freq)) + geom_histogram() + theme_classic() + ylab("Count") + xlab("Cells per segregant")
ggsave("figures/s1.png")


ap_match = readRDS("../rproj/out/combined/Ap/segMatchList.RDS")
ap_match_round =  ifelse(ap_match$vg> 0.5,1,0)
ap_match_round2 = ifelse(ap_match$sprevG > 0.5,1,0)

ap_match_combined = t(ap_match_round2[,ap_match$best_match_seg])
#ap_match-
#i#ap_match_combined

accuracy = apply(ap_match_combined == ap_match_round, 1, sum)/ncol(ap_match_combined)

ap_df = readRDS("../rproj/out/combined/Ap/segData.RDS")
ap_df$barcode.features$accuracy = accuracy

ap_df$barcode.features %>% ggplot(aes(y=accuracy * 100,x=(nUMI))) + geom_point() +  scale_x_log10() + xlab(expression('UMI count')) + ylab("Genotype Accuracy") + theme_bw() + stat_cor(method="spearman",size=14) +
  theme(text=element_text(size=24)) + ylim(c(50,100)) + geom_hline(yintercept = median(ap_df$barcode.features$accuracy*100),color="red")

ggsave("figures/s2.png")

p3 = ap_combined %>% mutate(single_cell_sig = p_adj.old < 0.05) %>% 
  ggplot(aes(x=beta.old,y=beta.hmm,color=factor(single_cell_sig))) + geom_point(size=1.5) + 
  theme_bw()  + stat_cor(method="spearman",size=8,cor.coef.name = "rho",show.legend = F) + theme(text=element_text(size=18))  + ylab("local eQTL effect (bulk)") + 
  xlab("local eQTL effect (one-pot)") + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  ylim(c(-3,3))  + xlim(c(-4,4))  +  scale_color_manual(name="One-pot FDR",labels=c(">= 0.05","< 0.05"),values=c("#e6ab02","#1b9e77")) +
  guides(color = guide_legend(override.aes = list(size = 8)))
p3
