pb = mating_efficiency_gpa1_strains%>% ggplot(aes(y=norm,x=strain_new)) + geom_boxplot(outlier.shape  = NA) +  geom_jitter(width=.2,size=3) +
  ylab("Mating efficiency") + xlab("Strain") + facet_wrap(~color2)+ theme_bw()  +
  theme(text=element_text(size=18))  +stat_compare_means(comparisons = mycomparsions,method="t.test",size=8) + theme(text=element_text(size=18))#
pb
ggsave("fig_final/s9.png",dpi=300,width=16,height=12)

