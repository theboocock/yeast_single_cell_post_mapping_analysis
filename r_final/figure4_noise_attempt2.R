se_em = combined_objects$noise$ASE_EMMEANS

#A = combined_objects$noise$ASE %>% filter(cross == "A")
#cor(A$estimate.disp,A$estimate.cond,use = "pairwise.complete.obs")
#plot(A$estimate.cond,A$estimate)
ase_em2 = ase_em %>% filter(theta > -6) %>% filter(emmean.mean > -7) #%>% ggplot(aes(y=theta,x=emmean.mean)) + geom_point()
## ##
sig_disp = combined_objects$noise$ASE %>% filter(p_adj_disp < 0.05)
sig_mean = combined_objects$noise$ASE %>% filter(p_adj_ase < 0.05)

sig_simple_df = sig_disp %>% dplyr::select(gene,cross,p_adj_disp,p.value.disp) %>% mutate(gene_cross=paste(gene,cross,sep="_"))
sig_simple_df2 = sig_mean %>% dplyr::select(gene,cross) %>% mutate(gene_cross=paste(gene,cross,sep="_"))


#sig_disp

df2 = ase_em2 %>% group_by(gene,cross) %>% summarise(mean_shift=emmean.mean[2]-emmean.mean[1], disp_shit=theta[2]-theta[1],se_disp=mean(c(SE.disp[1],SE.disp[2])),se_mean=mean(c(SE.mean[1],SE.mean[2])))
df2 = df2 %>% mutate(disp_pos = ifelse(disp_shit > 0, disp_shit, - disp_shit), mean_pos=ifelse(disp_shit > 0, mean_shift,-mean_shift))#
df2 = df2 %>% mutate(gene_cross=paste(gene, cross, sep="_"))

df2 = df2 %>% mutate(is_sig_disp = ifelse(gene_cross %in% sig_simple_df$gene_cross,T,F)) %>% 
  mutate(is_sig_ase = ifelse(gene_cross %in% sig_simple_df2$gene_cross,T,F))

library(robustbase)
library(sandwich)
library(lmtest)
library(modelr)
library(broom)
df3 = df2  %>% filter(is_sig_disp | is_sig_ase)
m1 = lmrob(disp_shit ~ mean_shift,data=df3)
pred_df = predict(m1,newdata=df3,interval="confidence")

pred_df3 = cbind(pred_df,df3)
pred_df3_raw = cbind(pred_df_raw,df3)
pred_df3 %>% ggplot(aes(y=fit,x=mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70") + geom_line()  + 
  geom_point(data=df3,aes(x=mean_shift,y=disp_shit)) #+ geom_line(aes(y=fit,x=mean_shift,color="red"),data=pred_df3_raw)

df3 = df3 %>% mutate(low_mean=mean_shift - 1.96 * se_mean, high_mean = mean_shift - 1.96 * se_mean, low_disp=disp_shit - 1.96 * se_disp,high_disp=disp_shit + 1.96 * se_disp)


pred_df3 %>% ggplot(aes(y=fit,x=mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70") + geom_line()  + 
  geom_line(aes(y=fit,x=mean_shift),data=pred_df3)

pred_df4 = cbind(pred_df,df3)


#outside_trend = 
pred_df4$good = rep(NA,nrow(pred_df4))
pred_df4$dist = rep(NA,nrow(pred_df4))
for(i in 1:nrow(pred_df4)){
  print(i)
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
  else if(disp_lwr >  upr){
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
  geom_line(aes(y=fit,x=mean_shift),data=pred_df3)

pred_df5 = pred_df5 %>% left_join(sig_simple_df,by=c("gene_cross")) %>% mutate(log10p=-log10(p.value.disp))
pred_df5 = pred_df5 %>% inner_join(gene_expr_1000_df,by=c("gene.x"="transcript"))
pred_df5 %>% ggplot(aes(y=fit,mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70")  +geom_line()+ 
  geom_point(data=pred_df5,aes(x=mean_shift,y=disp_shit,color=good_disp,size=log10p)) + theme_bw() +
  xlab("ln(expression) delta") + ylab("ln(dispersion) delta")  # + geom_line()  +

pred_df5 %>% ggplot(aes(y=fit,mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70")  +geom_line()+ 
  geom_point(data=pred_df5,aes(x=mean_shift,y=disp_shit,color=good_disp)) + theme_bw() +
  xlab("ln(expression) delta") + ylab("ln(dispersion) delta")  # + g

#install.packages("olsrr")
#library(olsrr)

(pred_df5)



pred_df4$resid = residuals(m1)
hist(residuals(m1),breaks=10)
plot(pred_df4$fit,pred_df4$disp_shit)
qqresid(m1)

pred_df4$idx = 1:nrow(pred_df4)
resid_df = data.frame(resid=(residuals(m1)),idx=as.numeric(names(residuals(m1))))

blah =pred_df4 %>% left_join(resid_df,by=("idx"))

blah %>% filter(abs(resid)> 2)

plot(residuals(m1))

#pred_df4

o#ls_plot_resid_qq(m1)

pred_df5 %>% filter(gene.x == "YFL014W")


pred_df4 = pred_df4 %>% mutate(label=ifelse(gene == "YMR303C",paste(gene,cross,sep="_"),NA))
pred_df3  %>% mutate(label=ifelse(gene == "YFL014W",paste(gene,cross,sep="_"),NA)) %>% ggplot(aes(y=fit,x=mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70") + geom_line()  + 
  geom_point(data=pred_df4,aes(x=mean_shift,y=disp_shit,color=good)) + geom_label_repel(aes(x=mean_shift, y=disp_shit, label=label),box.padding = 2,max.overlaps = Inf,size=6,data=pred_df4)# + geom_errorbar(data=pred_d