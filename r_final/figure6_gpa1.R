p_growth = growth_gpa_strains %>% mutate(plate_l = grepl("Les",plate)) %>% 
  ggplot(aes(y=resid,x=strain_new)) + geom_boxplot() + geom_jitter(width=0.2) + theme_bw() + ylab("Doublings per hour") + xlab("Strain")  + 
  theme(text=element_text(size=18)) 
p_mating = mating_efficiency_gpa1_strains%>% ggplot(aes(y=norm,x=strain_new)) + geom_boxplot(outlier.shape  = NA) +  geom_jitter(width=.2,size=3) +
  ylab("Mating efficiency") + xlab("Strain") + scale_color_manual(values=c("blue","red")) + theme_bw()  + # +
  theme(text=element_text(size=18)) # +stat_compare_means(comparisons = mycomparsions,method="t.test",size=8)#
plot_grid(p_growth,p_mating,nrow=2)
ggsave("fig_final/main/staging/figure6_left.svg",width=8,height=12)

#aa_2@phylo

#plot_grid(p_growth,p_mating,nrow=2)
aa_2 = as.treedata(mid)
gfive = groupOTU(aa_2@phylo, list(het_df$sample[het_df$gpa == 2],het_df$sample[het_df$gpa == 1],het_df$sample[het_df$gpa == 0]))
p3 = ggtree(gfive, aes(color=group), layout="fan",open.angle = 180)  + scale_color_manual(values=c("#377eb8","#e41a1c","grey85")) #

rotate_tree(p3,angle= 90) #+ scale_x_reverse()
ggsave("fig_final/main/staging/figure6_right_flipped_tree.svg",width=8,height=12)

#ggtree(groupOTU(aa_2@phylo,list(het_df$sample[het_df$Clades_trim == "3. Brazilian bioethanol"],het_df$sample[het_df$Clades_trim == "26. Asian fermentation"])),aes(color=group),layout="fan",open.angle = 180)


#ggtree(groupOTU(aa_2@phylo,list(het_df$sample[het_df$Clades_trim == "8. Mixed origin"],het_df$sample[het_df$Clades_trim == "3. Brazilian bioethanol"],het_df$sample[het_df$Clades_trim == "26. Asian fermentation"])),aes(color=group),layout="fan",open.angle = 180)




#list_group= list()
#i = 1
#for(mosaic in unique(het_df$Clades_trim[het_df$mosaic])){
#  print(mosaic)
#  list_group[[i]] = het_df$sample[het_df$Clades_trim == mosaic]
#  i = i + 1
#}
#
#list_group= list()
#i = 1
#for(mosaic in unique(het_df$Clades_trim)){
#  print(mosaic)
#  list_group[[i]] = het_df$sample[het_df$Clades_trim == mosaic]
#  i = i + 1
#}
#names(list_group) = unique(het_df$Clades_trim)


#rotate_tree(ggtree(groupOTU(aa_2@phylo,list_group),aes(color=group),layout="fan",open.angle = 180),angle=90)


#rotate_tree(ggtree(groupOTU(aa_2@phylo,list(het_df$sample[het_df$mosaic],het_df$sample[!het_df$mosaic])),aes(color=group),layout="fan",open.angle = 180), angle=90)

# #






