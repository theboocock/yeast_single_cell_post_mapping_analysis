library(tidyverse)
library(ggpubr)
library(cowplot)
library(viridis)
library(Seurat)
source("vars.R")
#source("load_cross_objects.R")

dir.create("figures")
cis_test_1000 = readRDS("../data/bulkEQTL_cis_only_test.RDS")
cis_test_1000$p_adj = p.adjust(cis_test_1000$p.value,method="fdr")
sum(cis_test_1000$p_adj < 0.05)
ap_cis = readRDS("../rproj/out/combined/Ap/geno_cis_comp.RDS")

ap_cis$prevGenoCis = ap_cis$prevGenoCis %>% mutate(p_adj=p.adjust(p,method = "fdr"))
ap_cis$hmmGenoCis = ap_cis$hmmGenoCis %>% mutate(p_adj=p.adjust(p,method = "fdr"))
ap_combined = ap_cis$prevGenoCis %>% inner_join(ap_cis$hmmGenoCis,by=c("gene"),suffix=c(".old",".hmm"))
ap_combined = ap_combined %>% mutate(sig_both = beta.old * beta.hmm > 0 & p_adj.old < 0.05 & p_adj.hmm < 0.05)
ap_combined_1000 = ap_combined %>% full_join(cis_test_1000,by=c("gene"="transcript"))
#ap_combined %>% ggplot(aes(y=beta.old,x=))
one_pot_bulk =cross_data$A$cis_test_with_disp_expr %>% full_join(cis_test_1000,by=c("transcript"="transcript"))

sum(one_pot_bulk$FDR < 0.05 & one_pot_bulk$Beta * one_pot_bulk$coefficient > 0,na.rm=T)


p1=  magick::image_read_svg("../figure_and_tables_paper/Figures/393.svg",width = 1600)
p2 = ggdraw() + draw_image(p1)
p3 = one_pot_bulk %>% mutate(single_cell_sig = FDR < 0.05) %>%filter(!is.na(Beta) & !is.na(coefficient)) %>% 
  ggplot(aes(x=Beta,y=coefficient,color=factor(single_cell_sig))) + geom_point(size=1.5) + 
  theme_bw()  + stat_cor(method="spearman",size=8,cor.coef.name = "rho",show.legend = F) + theme(text=element_text(size=18))  + ylab("local eQTL effect (bulk)") + 
  xlab("local eQTL effect (one-pot)") + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  ylim(c(-3,3))  + xlim(c(-4,4)) + 
  scale_color_brewer(name="One-pot FDR",labels=c(">= 0.05","< 0.05"),palette ="Dark2")  + 
  guides(color = guide_legend(override.aes = list(size = 8)))#guides(override.aes = aes(label = ""))
 # scale_color_discrete()
pcompare = ap_cis$hmmGenoCis %>% mutate(single_cell_sig = p_adj < 0.05) %>% 
  inner_join(ap_cis$prevGenoCis,by=c("gene")) %>% ggplot(aes(y=beta.y,x=beta.x,color=single_cell_sig)) + 
  geom_point(size=1.5) +
  geom_abline() + theme_bw()  + stat_cor(method="spearman",size=8,show.legend = F) +
  theme(text=element_text(size=18))  + xlab("local eQTL effect (single-cell HMM)") + 
  ylab("local eQTL effect (matched whole genome)") + #+ theme_ +
scale_color_brewer(name="HMM FDR",labels=c(">= 0.05","< 0.05"),palette ="Dark2") + # + 
  guides(color = guide_legend(override.aes = list(size = 8)))
library(patchwork)
#p2 + (pcompare/p3) + plot_annotation(tag_levels = "A")


p1=  magick::image_read_svg("../figure_and_tables_paper/Figures/one_pot.svg",width = 2000)
p4 = ggdraw() + draw_image(p1)

plot_grid(p4,nrow=2)

p_one_pot = cowplot::ggdraw() + cowplot::draw_image("../figure_and_tables_paper/Figures/one_pot.png")




by_new = cross_data$A$trans$hotspot_table
by_new$new_chrom = convert_chrom_to_simple_factor(by_new$chr)
by_old = cross_data$A_bulk$trans$hotsot_table
by_old$new_chrom = convert_chrom_to_simple_factor(by_old$chr)

p2 = ggplot(by_new,aes(pos,count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
  geom_bar(stat="identity", width=50000)+
  xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
  facet_grid(~new_chrom,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme_bw() +  theme(panel.spacing=unit(0, "lines"))

p3 = ggplot(by_old,aes(pos,count)) + #, col=LOD))+scale_colour_viridis_b()+ #LOD))) +#pha=-log10(LOD)/6))+geom_point()+
  geom_bar(stat="identity", width=50000)+
  xlab('')+ylab('')+ scale_alpha(guide = 'none') + 
  facet_grid(~new_chrom,scales="free", space="free")+        #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme_bw() +  theme(panel.spacing=unit(0, "lines"))

plot2 = p3 + geom_hline(yintercept = cross_data$A_bulk$trans$hotspot_threshold,color="red") +  theme_classic() + 
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +  theme(plot.margin = unit(c(0,0,0,0), "cm"))   +
  scale_y_reverse(expand=c(0,0)) + 
  xlab("Genome positiion (Mb)") + ylab("Distal eQTLs") +
  scale_x_continuous(breaks = c(.5e6,1e6,1.5e6,2e6),labels = function(x)x/1e6)  + 
  theme(text=element_text(size=18))  + ylab("Distal eQTLs (bulk)")  + 
  coord_cartesian(ylim=c(1500,0)) 
  
plot1 = p2 + geom_hline(yintercept = cross_data$A$trans$hotspot_threshold,color="red")   + theme_classic() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +  
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) + scale_y_continuous(expand=c(0,0)) +
  theme(text=element_text(size=18)) + coord_cartesian(ylim=c(0,100)) + ylab("Distal eQTLs (single-cell)")# + ggtitle("Single-cell eQTL")
  
#plot1 / plot2

pa = cross_data$A$trans$combined_peaks %>% filter(FDR < 0.05)  %>% 
  sample_frac() %>% filter(tchr!= "chrmt") %>%
  ggplot(aes(x=pos,y=tpos)) + 
  geom_point() + 
  facet_grid(tchrom_short_f~chrom_short_f,scales="free",space = "free",switch="both") + 
  theme_classic() + scale_color_brewer(palette ="Dark2") + theme(text=element_text(size=18)) + 
  ylab("Transcript position (Mb)") + xlab("Variant position (Mb)") + 
  scale_x_continuous(breaks = c(.5e6,1e6,1.5e6,2e6),labels = function(x)x/1e6) +
  scale_y_continuous(breaks =c(.5e6,1e6,1.5e6,2e6), labels =  function(x)x/1e6) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "grey20"))
  




#p4 / (pcompare+p3) + plot_annotation(tag_levels = "A") + plot_layout(height=c(1,2))
#ggsave("figures/figure1.png")


#pg = plot_grid(pcompare,p3,labels=c("B","C"),label_size = 18)
#cowplot::plot_grid(p4,pg,labels=c("A",""),label_size = 18,nrow=2,rel_heights = c(1,2))
#ggsave("figures/figure1.png",dpi=300)

bot_right = cowplot::plot_grid(plot1,plot2,nrow=2,labels=c("C","D"),label_size = 18,align = "v")
bot = cowplot::plot_grid(pa,bot_right,labels=c("B",""),label_size = 18)
cowplot::plot_grid(p_one_pot,bot,rel_heights = c(1,2.5),nrow=2,labels = c("A",""),label_size = 18)

#cowplot::plot_grid(p_one_pot,pa + (plot1 / plot2),nrow=2,rel_heights = c(1,3))
ggsave("figures/figure1.png",dpi=300,bg="white")

ggsave("figures/figure1.svg",dpi=300,bg="white")
