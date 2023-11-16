
#source("vars.R")
source("load_cross_objects.R")

#dir.create("figures")
cis_test_1000 = readRDS("data/bulkEQTL_cis_only_test.RDS")
cis_test_1000$p_adj = p.adjust(cis_test_1000$p.value,method="fdr")
sum(cis_test_1000$p_adj < 0.05)
ap_cis = readRDS("data/out/combined/Ap/geno_cis_comp.RDS")

ap_cis$prevGenoCis = ap_cis$prevGenoCis %>% mutate(p_adj=p.adjust(p,method = "fdr"))
ap_cis$hmmGenoCis = ap_cis$hmmGenoCis %>% mutate(p_adj=p.adjust(p,method = "fdr"))
ap_combined = ap_cis$prevGenoCis %>% inner_join(ap_cis$hmmGenoCis,by=c("gene"),suffix=c(".old",".hmm"))
ap_combined = ap_combined %>% mutate(sig_both = beta.old * beta.hmm > 0 & p_adj.old < 0.05 & p_adj.hmm < 0.05)
ap_combined_1000 = ap_combined %>% full_join(cis_test_1000,by=c("gene"="transcript"))
#ap_combined %>% ggplot(aes(y=beta.old,x=))
one_pot_bulk =cross_data$A$cis$cis_test_with_disp_expr %>% full_join(cis_test_1000,by=c("transcript"="transcript"))

sum(one_pot_bulk$FDR < 0.05 & one_pot_bulk$Beta * one_pot_bulk$coefficient > 0,na.rm=T)

nrow(ap_combined %>% filter(p_adj.old < 0.05))
nrow(ap_combined %>% filter(p_adj.hmm < 0.05))

#p1=  magick::image_read_svg("../figure_and_tables_paper/Figures/one_p",width = 1600)
#p1 = ggdraw() + draw_image(p1)


p1 = ggdraw() + draw_image("data/images//one_pot.png")

p3 = one_pot_bulk %>% mutate(single_cell_sig = FDR < 0.05) %>%filter(!is.na(Beta) & !is.na(coefficient)) %>% 
  ggplot(aes(x=Beta,y=coefficient,color=factor(single_cell_sig))) + geom_point(size=1.5) + 
  theme_bw()  + theme(text=element_text(size=18))  + ylab("local eQTL effect (bulk)") + 
  xlab("local eQTL effect (one-pot)") + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  ylim(c(-1,1))  + xlim(c(-1,1))  +  scale_color_manual(name="One-pot FDR",labels=c(">= 0.05","< 0.05"),values=c("#e6ab02","#1b9e77")) +
  guides(color = guide_legend(override.aes = list(size = 8)))#guides(override.aes = aes(label = ""))
rds = readRDS("../rproj/out/cell_cycle/00_BYxRM_480MatA_1/cell_cycle_final.RDS")
p2 = DimPlot(rds$cc_seurat, group.by="cell_cycle",label = T,label.size = 12,pt.size = 2) + theme(text=element_text(size=18)) + ggtitle("") + scale_color_brewer(palette = "Set2") + theme(legend.position = "none")
pg = plot_grid(p2,p3,labels=c("B","C"),label_size = 18)
plot_grid(p1,pg,labels=c("A",""),nrow=2,rel_heights = c(1,2))
ggsave("fig_final//figure1.png",bg="white",dpi=300,width=16,height=12)
ggsave("fig_final//figure1.svg",bg="white",dpi=300,width=16,height=12)

#### Figure 1, done #####




