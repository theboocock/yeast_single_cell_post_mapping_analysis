seg_match = readRDS("data/out/combined/Ap/segMatch.RDS")

seg_match_df = data.frame(seg=names(table(seg_match$best_match_seg)),table(seg_match$best_match_seg))
med_seg = median(seg_match_df$Freq)
seg_match_df %>% ggplot(aes(x=Freq)) + geom_histogram() + theme_classic() + ylab("Count") + xlab("Cells per segregant") + theme(text=element_text(size=24)) + geom_vline(xintercept = med_seg,color="red")
ggsave("fig_final//s6.png",width=16,height=12,bg="white",dpi=300)
ggsave("fig_final/svg/s6.svg",width=16,height=12,bg="white",dpi=300)
