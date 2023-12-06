ase_em = combined_objects$noise$ASE_EMMEANS

#A = combined_objects$noise$ASE %>% filter(cross == "A")
#cor(A$estimate.disp,A$estimate.cond,use = "pairwise.complete.obs")
#plot(A$estimate.cond,A$estimate)
ase_em2 = ase_em #%>% filter(theta > -6) %>% filter(emmean.mean > -7) #%>% ggplot(aes(y=theta,x=emmean.mean)) + geom_point()


## ##
sig_disp = combined_objects$noise$ASE %>% filter(p_adj_disp < 0.05)
sig_mean = combined_objects$noise$ASE %>% filter(p_adj_ase < 0.05)
combined_objects$noise$ASE = combined_objects$noise$ASE  %>% mutate(gene_cross=paste(gene,cross,sep="_"))
sig_simple_df = sig_disp %>% dplyr::select(gene,cross,estimate.cond,estimate.disp,p.value.cond,p.value.disp,p_adj_ase, p_adj_disp) %>% mutate(gene_cross=paste(gene,cross,sep="_"))
sig_simple_df2 = sig_mean %>% dplyr::select(gene,cross) %>% mutate(gene_cross=paste(gene,cross,sep="_"))
#sig_disp

df2 = ase_em2 %>% group_by(gene,cross) %>% summarise(mean_shift=emmean.mean[2]-emmean.mean[1], disp_shit=theta[2]-theta[1],se_disp=mean(c(SE.disp[1],SE.disp[2])),se_mean=mean(c(SE.mean[1],SE.mean[2])))
df2 = df2 %>% mutate(disp_pos = ifelse(disp_shit > 0, disp_shit, - disp_shit), mean_pos=ifelse(disp_shit > 0, mean_shift,-mean_shift))#
df2 = df2 %>% mutate(gene_cross=paste(gene, cross, sep="_"))

df2 = df2 %>% mutate(is_sig_disp = ifelse(gene_cross %in% sig_simple_df$gene_cross,T,F)) %>% 
  mutate(is_sig_ase = ifelse(gene_cross %in% sig_simple_df2$gene_cross,T,F))

df3 = df2  %>% filter(is_sig_disp | is_sig_ase)
df3 = df3 %>% left_join(combined_objects$noise$ASE,by=c("gene_cross"))
m1 = lmrob(theta ~ estimate.cond,data=df3)
pred_df = predict(m1,newdata=df3,interval="confidence")
pred_df3 = cbind(pred_df,df3)
df3 = df3 %>% mutate(low_mean=estimate.cond - 1.96 * std.error.cond, 
                     high_mean = estimate.cond + 1.96 * std.error.cond, 
                     low_disp=theta - 1.96 * std.error.disp,
                     high_disp=theta + 1.96 * std.error.disp)


#df3 %>% ggplot(aes(y=disp_shit,x=estimate.disp)) + geom_point()

pred_df4 = cbind(pred_df,df3)


#outside_trend = 
pred_df4$good = rep(NA,nrow(pred_df4))
pred_df4$dist = rep(NA,nrow(pred_df4))
for(i in 1:nrow(pred_df4)){
  #print(i)
  lwr = pred_df4$lwr[i]
  upr = pred_df4$upr[i]
  disp_lwr = pred_df4$low_disp[i]
  disp_upr = pred_df4$high_disp[i]
  
  #disp_lwr
  good = F
  if(is.na(disp_lwr) | is.na(disp_upr)){
    good = F
    dist = 0
  }
  else if(disp_lwr > upr){
    good = T
    dist = abs(disp_lwr - upr)
  }else if(disp_upr < lwr){
    good = T
    dist = abs(disp_upr - lwr)
  }
  pred_df4$good[i] = good
  pred_df4$dist[i] = dist
  
}
pred_df5= pred_df4 %>% mutate(good_disp=good & is_sig_disp) 



#pred_df3 %>% ggplot(aes(y=fit,x=mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70") + geom_line()  + 
#  geom_line(aes(y=fit,x=mean_shift),data=pred_df3)

#pred_df5 = pred_df5 %>% left_join(combined_objects$noise$ASE,by=c("gene_cross")) %>% mutate(log10p=-log10(p.value.disp))
pred_df5 = pred_df5 %>% left_join(gene_expr_1000_df,by=c("gene.x"="transcript"))
#pred_df5 %>% ggplot(aes(y=fit,mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70")  +geom_line()+ 
#  geom_point(data=pred_df5,aes(x=mean_shift,y=disp_shit,color=good_disp,size=log10p)) + theme_bw() +
#  xlab("ln(expression) delta") + ylab("ln(dispersion) delta")  # + geom_line()  +
#c(-4)
pred_df5 = pred_df5 %>% inner_join(genes_name_trans,by=c("gene.x"="gene_id"))
#pred_df5 = pred_df5 %>% mutate(hsp_label=ifelse(gene_name=="HSP12",paste("HSP12"," ",cross.x),""))

p2 = pred_df5 %>% ggplot(aes(y=fit,x=estimate.cond)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70")  +geom_line()+ 
  geom_point(data=pred_df5,aes(x=estimate.cond,y=theta,color=good_disp)) + theme_bw() +
  xlab("ln(expression) delta") + ylab("ln(dispersion) delta")  + xlim(c(-5,5) ) + ylim(c(-5,5))  + theme(text=element_text(size=18)) +
  scale_color_manual(name = "95% noise confidence interval\n of noise\noverlaps trend line",labels=c("Yes","No"),values=c("grey60","red")) + ylab(expression(Delta*ln*"(noise)")) + 
  xlab(expression(Delta*ln*"(expression)")) + theme(legend.position = "none") #+ geom_text() #+ geom_hline(yintercept = 0) + geom_vline(xintercept = 0) 

#+ ggrepel::geom_label_repel(aes(x=mean_shift,y=disp_shit,label=hsp_label),max.overlaps = Inf,size=12,min.segment.length = 3)  
  #theme(title = "min.segment.length = 0"


#p2 = pred_df5 %>% ggplot(aes(y=fit,mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70")  +geom_line()+ 
#  geom_point(data=pred_df5,aes(x=mean_shift,y=disp_shit,color=good_disp)) + theme_bw() +
#  xlab("ln(expression) delta") + ylab("ln(dispersion) delta")  + theme(text=element_text(size=18)) +
#  scale_color_manual(name = "95% noise confidence interval\n of noise\noverlaps trend line",labels=c("Yes","No"),values=c("grey60","red")) + ylab(expression(Delta*ln*"(noise)")) + 
#  xlab(expression(Delta*ln*"(expression)")) + theme(legend.position = "none") #+ geom_text() #+ geom_hline(yintercept = 0) + geom_vline(xintercept = 0) 


#p2 = pred_df5 %>% ggplot(aes(y=fit,mean_shift-fit)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70")  +geom_line()+ 
#  geom_point(data=pred_df5,aes(x=mean_shift-fit,y=disp_shit,color=good_disp)) + theme_bw() +
#  xlab("ln(expression) delta") + ylab("ln(dispersion) delta")  + xlim(c(-2,2) ) + ylim(c(-2,2))  + theme(text=element_text(size=24)) +
#  scale_color_discrete(name = "Overlaps the\ntrend line",labels=c("No","Yes")) + ylab(expression(Delta*ln*"(noise)")) + 
#  xlab(expression(Delta*ln*"(expression)")) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  

tc = rowSums(cross_data$A$cis$ASE$ase_data$counts)/ncol(cross_data$A$cis$ASE$ase_data$counts)
set.seed(1)
mu_big = quantile(tc,prob=.99)
n1 = rnegbin(n=5000,mu=mu_big,theta = 1000)
n2 = rnegbin(n=5000,mu=mu_big,theta = 1)
d1 = data.frame(count=n1,allele="A")
d2 = data.frame(count=n2,allele="B")
df_theta_diff = rbind(d1,d2) 

p_theta = df_theta_diff  %>% ggplot(aes(x=count,color=allele)) + stat_ecdf() + theme_bw()  + theme(text=element_text(size=18)) + ylab("Probability") + 
  scale_color_brewer(name="Allele",palette = "Dark2")  + coord_cartesian(xlim=c(0,40))#+ facet_wrap(~allele)

#p14 = df %>% ggplot(aes(x=count,color=allele)) + stat_ecdf()#+ facet_wrap(~allele)
shift = quantile(exp(pred_df5$mean_shift),na.rm=T,prob=.9)
n1 = rnegbin(n=5000,mu=mu_big,theta = 1)
n2 = rnegbin(n=5000,mu=mu_big+shift,theta = 1)

d1 = data.frame(count=n1,allele="A")
d2 = data.frame(count=n2 ,allele="B")
df_mean = rbind(d1,d2) 
p_mean = df_mean  %>% ggplot(aes(x=count,color=allele)) + stat_ecdf() + theme_bw()  + theme(text=element_text(size=18)) + ylab("Probability") + 
  scale_color_brewer(name="Allele",palette = "Dark2") + coord_cartesian(xlim=c(0,40))  #+ facet_wrap(~allele)

pa  = plot_grid(p_mean,p_theta,ncol=2,labels=c("A","B"))


plot_grid(pa,p2,nrow=2,rel_heights = c(1,3),labels=c("","C"),label_size = 12)

ggsave("fig_final/main/staging///figure4.png",dpi=300,width=16,height=12,device = png)
ggsave("fig_final/main/staging/figure4.svg",dpi=300,width=16,height=12)

#p14 = df %>% ggplot(aes(x=count,color=allele)) + stat_ecdf()#+ facet_wrap(~allele)
#m3 = glmmTMB::glmmTMB(df$count ~ df$allele, dispformula = ~df$allele,family= glmmTMB::nbinom2())
#p15 = df %>% ggplot(aes(x=count,color=allele)) + geom_density()#+ facet_wrap(~allele)

#p14
#p15

# + g

#install.packages("olsrr")
#library(olsrr)

#