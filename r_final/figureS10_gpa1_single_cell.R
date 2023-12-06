p1 = DimPlot(gpa_c, group.by="cell_cycle",label=T,label.size = 18) + ggtitle("")
# A  =420/421 82R
# B = 416/417 WT
gpa1=gpa_c@meta.data$dataset

gpa_c@meta.data$cell_cycle2 = convert_cell_cycle_to_simple_factor(gpa_c@meta.data$cell_cycle)

gpa1[gpa1=="A"] = "82R"
gpa1[gpa1=="B"] = "WT"
gpa_c@meta.data$gpa1 = gpa1
p2 = DimPlot(gpa_c, group.by="gpa1",pt.size = .5) + ggtitle("")
#gpa_c@meta.data
#table(gpa_c$cell_cycle == "G1",gpa_c$gpa1)  %>% as_data_frame()%>% pivot_longer()
#plot_grid(p1,p2)
df_gpa1 = gpa_c@meta.data %>% group_by(cell_cycle2) %>% summarise(gpa1=sum(gpa1=="82R"),n=n(),freq=gpa1/n())
df_gpa1$se = sqrt(df_gpa1$freq*(1-df_gpa1$freq)/df_gpa1$n)

p3 = df_gpa1 %>% ggplot(aes(y=freq,x=cell_cycle2))+ geom_point(stat="identity") + geom_errorbar(aes(ymax=freq + 1.96*se, ymin=freq - 1.96*se))+ ylab("82R Frequency") + xlab("cell-cycle stage") +
  theme_bw() + theme(text=element_text(size=24))

plot_grid(p1,p2,p3,ncol=3,label_size = 24, labels=c("A","B","C"))
ggsave("fig_final/s10.png",width=16,dpi=300,height=12)
ggsave("fig_final/svg/s10.svg",width=16,dpi=300,height=12)
