het_df$Clades_trim[het_df$Clades_trim == ""]  = "Unassigned"
het_df %>% group_by(Clades_trim) %>% summarise(gpac=sum(gpa)/(n()*2), n=n())%>% mutate(clade_trime2 = paste(Clades_trim," (",n,")",sep=""))  %>%
  ggplot(aes(x=reorder(clade_trime2,gpac),y=gpac)) + geom_bar(stat="identity")   + coord_flip() + theme_classic() +
  theme(text=element_text(size=18)) + xlab("Clade") + ylab("82R allele-frequency")
ggsave(width=16,height=12,dpi=300,filename = "fig_final/s14.png")
ggsave(width=16,height=12,dpi=300,filename = "fig_final/svg/s14.svg")
