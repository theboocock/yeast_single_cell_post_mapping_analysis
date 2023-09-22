het_df_filt %>% ggplot(aes(y=Proportion.of.clean.heterozygous.SNPs..whole.dataset.,x=factor(gpa))) + geom_violin()+ geom_jitter(width=0.2)  + theme_classic() + 
 ylab("Proportion of sites called as heterozygous") + xlab("82R genotype") + theme(text=element_text(size=18))
ggsave("fig_final/s11.png",dpi=300,width=16,height=12)
