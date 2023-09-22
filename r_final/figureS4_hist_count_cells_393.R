seg_match = readRDS("../rproj/out/combined/Ap/segMatch.RDS")
seg_match_df = data.frame(seg=names(table(seg_match$best_match_seg)),table(seg_match$best_match_seg))
seg_match_df %>% ggplot(aes(x=Freq)) + geom_histogram() + theme_classic() + ylab("Count") + xlab("Cells per segregant")
ggsave("fig_final//s4.png",width=16,height=16,bg="white",dpi=300)
