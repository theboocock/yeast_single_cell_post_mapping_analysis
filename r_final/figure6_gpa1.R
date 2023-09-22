mycomparsions = list(c("WT","82R"),c("82R","469I"),c("WT","469I"))
p_growth = growth_gpa_strains %>% mutate(plate_l = grepl("Les",plate)) %>% 
  ggplot(aes(y=resid,x=strain_new)) + geom_boxplot() + geom_jitter(width=0.2) + theme_bw() + ylab("Doublings per hour") + xlab("Strain")  + theme(text=element_text(size=18)) +
  stat_compare_means(comparisons = mycomparsions,method="t.test",size=8) 
p_mating = mating_efficiency_gpa1_strains%>% ggplot(aes(y=norm,x=strain_new)) + geom_boxplot(outlier.shape  = NA) +  geom_jitter(width=.2,size=3) +
  ylab("Mating efficiency") + xlab("Strain") + scale_color_manual(values=c("blue","red")) + theme_bw()  +
  theme(text=element_text(size=18))  +stat_compare_means(comparisons = mycomparsions,method="t.test",size=8)#


