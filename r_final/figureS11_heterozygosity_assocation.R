het_df_filt = het_df_filt %>% mutate(gpa1_new = case_when(
                          gpa == 0 ~ "Trp/Trp",
                          gpa == 2 ~ "Arg/Arg"
                            ))
het_df_filt %>% ggplot(aes(y=Proportion.of.clean.heterozygous.SNPs..whole.dataset.,x=gpa1_new)) + geom_violin()+ geom_jitter(width=0.2)  + theme_classic() + 
 ylab("Proportion of sites called as heterozygous") + xlab("82R genotype") + theme(text=element_text(size=24)) + scale_x_discrete(limits = rev)
#+ scale_x_reverse()
ggsave("fig_final/s11.png",dpi=300,width=16,height=12)
