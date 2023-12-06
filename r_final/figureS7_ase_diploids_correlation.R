combined_objects$cis_eqtl_ase_only_with_onepot %>% filter(!is.na(Beta) & !is.na(estimate)) %>% ggplot(aes(x=Beta,y=estimate,color=FDR < 0.05)) + geom_point() + facet_wrap(~cross2,nrow=3) + theme_bw() + theme(text=element_text(size=24)) +
  coord_cartesian(ylim=c(-2,2)) +
  #stat_cor(method = "spearman") +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + xlab("local eQTL effect (one-pot)")  + ylab("allele-specific expression effect (diploid hybrids)") + 
  scale_color_brewer(name="One-pot FDR",labels=c(">= 0.05","< 0.05"),palette ="Dark2") +   
  guides(color = guide_legend(override.aes = list(size = 8)))


A = combined_objects$cis_eqtl_ase_only_with_onepot %>% dplyr::filter(cross == "A" & FDR < 0.05)
cor(A$estimate,A$Beta,method="spearman")
B = combined_objects$cis_eqtl_ase_only_with_onepot %>% dplyr::filter(cross == "B" & FDR < 0.05)
cor(B$estimate,B$Beta,method="spearman")
C = combined_objects$cis_eqtl_ase_only_with_onepot %>% dplyr::filter(cross == "3004" & FDR < 0.05)
cor(C$estimate,C$Beta,method="spearman")
#combined_objects$cis_eqtl_ase_noise %>% 


ggsave("fig_final//s7.png",bg="white", dpi=300, height=12,width=16)
ggsave("fig_final/svg//s7.svg", bg="white", dpi=300, height=12, width=16)

