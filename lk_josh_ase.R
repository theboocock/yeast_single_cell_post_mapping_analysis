
ase_em = combined_objects$noise$ASE_EMMEANS

#A = combined_objects$noise$ASE %>% filter(cross == "A")
#cor(A$estimate.disp,A$estimate.cond,use = "pairwise.complete.obs")
#plot(A$estimate.cond,A$estimate)
ase_em2 = ase_em %>% filter(theta > -6) %>% filter(emmean.mean > -7) #%>% ggplot(aes(y=theta,x=emmean.mean)) + geom_point()
## ##
sig_disp = combined_objects$noise$ASE %>% filter(p_adj_disp < 0.05)
sig_mean = combined_objects$noise$ASE %>% filter(p_adj_ase < 0.05)

sig_simple_df = sig_disp %>% dplyr::select(gene,cross) %>% mutate(gene_cross=paste(gene,cross,sep="_"))
sig_simple_df2 = sig_mean %>% dplyr::select(gene,cross) %>% mutate(gene_cross=paste(gene,cross,sep="_"))





m1 = (lm(theta ~ emmean.mean + cross,data=ase_em2))



df2 = ase_em2 %>% group_by(gene,cross) %>% summarise(mean_shift=emmean.mean[2]-emmean.mean[1], disp_shit=theta[2]-theta[1],se_disp=mean(c(SE.disp[1],SE.disp[2])),se_mean=mean(c(SE.mean[1],SE.mean[2])))
df2 = df2 %>% mutate(disp_pos = ifelse(disp_shit > 0, disp_shit, - disp_shit), mean_pos=ifelse(disp_shit > 0, mean_shift,-mean_shift))
df2 = df2 %>% mutate(gene_cross=paste(gene, cross, sep="_"))

df2 = df2 %>% mutate(is_sig_disp = ifelse(gene_cross %in% sig_simple_df$gene_cross,T,F)) %>% 
  mutate(is_sig_ase = ifelse(gene_cross %in% sig_simple_df2$gene_cross,T,F))

df2 %>% ggplot(aes(x=mean_shift,y=disp_shit,color=is_sig)) + geom_point(alpha=0.5) + facet_wrap(~cross) + scale_color_manual(values=c("grey50","red"))# + scale_color_brewer(alpha=0.5)


p1 = df2  %>% filter(is_sig_ase) %>% ggplot(aes(x=mean_shift,y=disp_shit)) + geom_point()  + ylab("ln(dispersion) delta") + xlab("ln(expression) delta") #+  stat_smooth(level=.99,method="lm") + 
 # geom_errorbar(aes(ymin=disp_shit- 1.96  * se_disp, ymax = disp_shit + 1.96 * se_disp))#$+ scale_color_manual(values=c("grey50","red"))# + scale_color_brewer(alpha=0.5)

p2 = df2  %>% filter(is_sig_disp) %>% ggplot(aes(x=mean_shift,y=disp_shit)) + geom_point()  + ylab("ln(dispersion) delta") + xlab("ln(expression) delta")

p3= df2  %>% filter(is_sig_disp | is_sig_ase) %>% ggplot(aes(x=mean_shift,y=disp_shit)) + geom_point()  + ylab("ln(dispersion) delta") + xlab("ln(expression) delta")


df2  %>% filter(!is_sig_disp & is_sig_ase) %>% ggplot(aes(x=mean_shift,y=disp_shit)) + geom_point()  + ylab("ln(dispersion) delta") + xlab("ln(expression) delta")


fisher.test(table(df2$is_sig_disp,df2$is_sig_ase))

cowplot::plot_grid(p1,p2,p3,ncol=3)

#install.packages("sandwich")
library(robustbase)
library(sandwich)
library(lmtest)
library(modelr)
library(broom)

df3 = df2  %>% filter(is_sig_disp)
m1 = lm(disp_shit ~ mean_shift ,data=df3)

pred_df_raw = predict(m1,newdata=df3,interval="confidence")

bptest(m1)

m1 = lmrob(disp_shit ~ mean_shift,data=df3)

pred_df = predict(m1,newdata=df3,interval="confidence")

pred_df3 = cbind(pred_df,df3)
pred_df3_raw = cbind(pred_df_raw,df3)
pred_df3 %>% ggplot(aes(y=fit,x=mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70") + geom_line()  + 
  geom_point(data=df3,aes(x=mean_shift,y=disp_shit)) + geom_line(aes(y=fit,x=mean_shift,color="red"),data=pred_df3_raw)

df3 = df3 %>% mutate(low_mean=mean_shift - 1.96 * se_mean, high_mean = mean_shift - 1.96 * se_mean, low_disp=disp_shit - 1.96 * se_disp,high_disp=disp_shit + 1.96 * se_disp)


pred_df4 = cbind(pred_df,df3)


#outside_trend = 
pred_df4$good = rep(NA,nrow(pred_df4))
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
  }
  else if(disp_lwr >  upr){
    good = T
  }else if(disp_upr < lwr){
    good = T
  }
  pred_df4$good[i] = good
  
}
#install.packages("olsrr")
library(olsrr)

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


pred_df4 = pred_df4 %>% mutate(label=ifelse(gene == "YMR303C",paste(gene,cross,sep="_"),NA))
pred_df3  %>% mutate(label=ifelse(gene == "YFL014W",paste(gene,cross,sep="_"),NA)) %>% ggplot(aes(y=fit,x=mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70") + geom_line()  + 
  geom_point(data=pred_df4,aes(x=mean_shift,y=disp_shit,color=good)) + geom_label_repel(aes(x=mean_shift, y=disp_shit, label=label),box.padding = 2,max.overlaps = Inf,size=6,data=pred_df4)# + geom_errorbar(data=pred_df4,aes(ymin=disp_shit- 1.96  * se_disp, ymax = disp_shit + 1.96 * se_disp)) 

pred_df3



names(table(pred_df4$gene[pred_df4$good]))[(table(pred_df4$gene[pred_df4$good]) > 2)]


pred_df3 %>% ggplot(aes(y=fit,x=mean_shift)) + geom_ribbon(aes(ymin=lwr,ymax=upr),fill="grey70") + geom_line()  + 
  geom_point(data=pred_df4,aes(x=mean_shift,y=disp_shit,color=good)) + geom_errorbar(data=pred_df4,aes(ymin=disp_shit- 1.96  * se_disp, ymax = disp_shit + 1.96 * se_disp)) #+ geom_line(aes(y=fit,x=mean_shift,color="red"),data=pred_df3_raw)


pred_df3 %>% mutate(label=ifelse(gene == "YFL014W",paste(gene,cross,sep="_"),NA))%>% ggplot(aes(y=disp_shit/se_disp,x=mean_shift/se_mean)) + geom_point() +  geom_label_repel(aes(label=label)) + stat_smooth(method="lm")

pred_df3 %>% filter()
summary(m1)

pred_df3 %>% filter(gene == "YFL014W")


library(ggrepel)
df2 %>% mutate(label=ifelse(gene == "YFL014W","HSP12",NA)) %>% filter(is_sig) %>% ggplot(aes(x=mean_shift,y=disp_shit)) + geom_point() + facet_wrap(~cross)  + ylab("ln(dispersion)") + xlab("ln(expression)") +  stat_smooth(level=.99,method="lm") + 
  geom_errorbar(aes(ymin=disp_shit- 1.96  * se_disp, ymax = disp_shit + 1.96 * se_disp)) + geom_errorbar(aes(xmin=mean_shift- 1.96*se_mean,xmax=mean_shift + 1.96*se_mean)) + geom_text_repel(aes(label=label))#$+ scale_color_manual(values=c("grey50","red"))# + scale_color_brewer(alpha=0.5)


df2_sig = df2  %>% filter(is_sig)

A = df2_sig %>% filter(cross == "A")
A = A %>% filter(!is.na(disp_shit)) %>% filter(!is.na(mean_shift))

am1  = lm(disp_shit ~ mean_shift, data=A)
#A = A %>% filter(!is.na(disp_shit))

amr = residuals(lm(disp_shit ~ mean_shift, data=A))

plot(amr, A$disp_shit)

am2 = predict(am,newdata=A)

plot(A$disp_shit,am2)
abline(0,1)
#df2 %>% left_join(sig_simple_df,by=c("gene","cross"))

df2 %>% filter(gene == "YFL014W")


plot(df2$disp_pos,df2$mean_pos)
plot(df2$disp_shit,df2$mean_shift)


df   = data.frame(disp_pos=df2$disp_pos, cross=df2$cross,mean_pos=df2$mean_pos)

m1 = lm(disp_pos ~ mean_pos+ cross,data=df2)





new_thetas= predict(m1,newdata=df)
plot(df2$disp_pos,new_thetas)
abline(0,1)

df$new_thetas = new_thetas

plot(df$new_thetas,df$theta_shit)

plot(df2$mean_shift,df2$disp_shit)

abline(0,1)


df
df



n1 = rnegbin(n=10000,mu=20,theta = 100000)
n2 = rnegbin(n=10000,mu=30,theta = 100000)

d1 = data.frame(count=n1,allele="A")
d2 = data.frame(count=n2,allele="B")
p2 = rbind(d1,d2) %>% ggplot(aes(x=count,color=allele)) + geom_density()#+ facet_wrap(~allele)
p1 = rbind(d1,d2) %>% ggplot(aes(x=count,color=allele)) + stat_ecdf()#+ facet_wrap(~allele)



n1 = rnegbin(n=100000,mu=20,theta = 7)
n2 = rnegbin(n=100000,mu=20,theta = 10000)

d1 = data.frame(count=n1,allele="A")
d2 = data.frame(count=n2,allele="B")
df = rbind(d1,d2) 
p3 = df %>% ggplot(aes(x=count,color=allele)) + stat_ecdf()#+ facet_wrap(~allele)
m3 = glmmTMB::glmmTMB(df$count ~ df$allele, dispformula = ~df$allele,family= glmmTMB::nbinom2())
p4 = df %>% ggplot(aes(x=count,color=allele)) + geom_density()#+ facet_wrap(~allele)



#df %>% ggplot(aes(x=))
#summary(m3)
#1/.69


n1 = rnegbin(n=100000,mu=20,theta = 7)
n2 = rnegbin(n=100000,mu=30,theta = 10000)

d1 = data.frame(count=n1,allele="A")
d2 = data.frame(count=n2,allele="B")
df = rbind(d1,d2) 
p5 = df %>% ggplot(aes(x=count,color=allele)) + stat_ecdf()#+ facet_wrap(~allele)
#m3 = glmmTMB::glmmTMB(df$count ~ df$allele, dispformula = ~df$allele,family= glmmTMB::nbinom2())
p6 = df %>% ggplot(aes(x=count,color=allele)) + geom_density()



n1 = rnegbin(n=100000,mu=20,theta = 10000)
n2 = rnegbin(n=100000,mu=30,theta = 7)

d1 = data.frame(count=n1,allele="A")
d2 = data.frame(count=n2,allele="B")
df = rbind(d1,d2) 
p7 = df %>% ggplot(aes(x=count,color=allele)) + stat_ecdf()#+ facet_wrap(~allele)
#m3 = glmmTMB::glmmTMB(df$count ~ df$allele, dispformula = ~df$allele,family= glmmTMB::nbinom2())
p8 = df %>% ggplot(aes(x=count,color=allele)) + geom_density()

#summary(m3)
cowplot::plot_grid(p1+ ggtitle("Diff mean same dispersion"),p2 + ggtitle("Diff mean same dispersion"),p3 + ggtitle("Diff dispersion, same mean"),p4 + ggtitle("Diff dispersion, same mean"),
                   p5 + ggtitle("higher mean, lower dispersion"),p6 + ggtitle("higher mean, lower dispersion"),p7 + ggtitle("higher mean, higher dispersion"),p8 + ggtitle("higher mean, higher dispersion"),ncol=4)


n1 = rnegbin(n=100000,mu=1,theta = 1000)
n2 = rnegbin(n=100000,mu=1,theta = 1)

d1 = data.frame(count=n1,allele="A")
d2 = data.frame(count=n2,allele="B")
df = rbind(d1,d2) 
p10 = df %>% ggplot(aes(x=count,color=allele)) + stat_ecdf()#+ facet_wrap(~allele)
m3 = glmmTMB::glmmTMB(df$count ~ df$allele, dispformula = ~df$allele,family= glmmTMB::nbinom2())
p11 = df %>% ggplot(aes(x=count,color=allele)) + geom_density()#+ facet_wrap(~allele)


n1 = rnegbin(n=100000,mu=1,theta = 1)
n2 = rnegbin(n=100000,mu=1.5,theta = 1)

d1 = data.frame(count=n1,allele="A")
d2 = data.frame(count=n2,allele="B")
df = rbind(d1,d2) 
p12 = df %>% ggplot(aes(x=count,color=allele)) + stat_ecdf()#+ facet_wrap(~allele)
m3 = glmmTMB::glmmTMB(df$count ~ df$allele, dispformula = ~df$allele,family= glmmTMB::nbinom2())
p13 = df %>% ggplot(aes(x=count,color=allele)) + geom_density()#+ facet_wrap(~allele)

plot_grid(p10,p11,p12,p13)


n1 = rnegbin(n=100000,mu=.1,theta = 1000)
n2 = rnegbin(n=100000,mu=.1,theta = 1)

d1 = data.frame(count=n1,allele="A")
d2 = data.frame(count=n2,allele="B")
df = rbind(d1,d2) 
p14 = df %>% ggplot(aes(x=count,color=allele)) + stat_ecdf()#+ facet_wrap(~allele)
m3 = glmmTMB::glmmTMB(df$count ~ df$allele, dispformula = ~df$allele,family= glmmTMB::nbinom2())
p15 = df %>% ggplot(aes(x=count,color=allele)) + geom_density()#+ facet_wrap(~allele)


n1 = rnegbin(n=100000,mu=1/10,theta = 1)
n2 = rnegbin(n=100000,mu=1.5/10,theta = 1)/

d1 = data.frame(count=n1,allele="A")
d2 = data.frame(count=n2,allele="B")
df = rbind(d1,d2) 
p16 = df %>% ggplot(aes(x=count,color=allele)) + stat_ecdf()#+ facet_wrap(~allele)
m3 = glmmTMB::glmmTMB(df$count ~ df$allele, dispformula = ~df$allele,family= glmmTMB::nbinom2())
p17 = df %>% ggplot(aes(x=count,color=allele)) + geom_density()#+ facet_wrap(~allele)


plot_grid(p14,p15,p16,p17)


aaa = rowSums(cross_data$`3004`$cis$ASE$ase_data$aC)
aaa[names(aaa) == "TDH3"]




pc = (cross_data$`3004`$cis$ASE$phased_counts[[1]] + cross_data$`3004`$cis$ASE$phased_counts[[2]])
dim(pc)
ax = colSums(pc)/nrow(pc)
ax[(ax> 30)]
hist((cross_data$`3004`$cis$ASE$geno_mean_disp$emmean.disp),breaks=1000)

cross_data$`3004`$cis$ASE$ase_noise %>% arrange(p_adj_disp)


cross_data$`3004`$cis$ASE$geno_mean_disp  %>% filter(gene == "YOR167C")


a = rnegbin(n=5000,mu=exp(2.969047),theta=exp(0.865195))
b = rnegbin(n=5000,mu=exp(2.861342),theta=exp(1.134289))
ag = data.frame(counts=c(a,b),allele=c(rep("A",length(a)),rep("B",length(b))))

ag %>% ggplot(aes(x=counts,color=allele)) + stat_ecdf()


norm = rowSums(cross_data$`3004`$cis$ASE$phased_counts[[1]]+ cross_data$`3004`$cis$ASE$phased_counts[[2]])

pc = cross_data$`3004`$cis$ASE$phased_counts[[1]][,"YOR167C"]# + cross_data$`3004`$cis$ASE$phased_counts[[2]])
pc2 = cross_data$`3004`$cis$ASE$phased_counts[[2]][,"YOR167C"]# + cross_data$`3004`$cis$ASE$phased_counts[[2]])


pcl = log1p(pc/(norm * 10000))
pcl2 = log1p(pc2/(norm * 10000))

df = data.frame(counts = c(pcl,pcl2),allele=c(rep("A",length(pc)),rep("B",length(pc2))))



df %>% ggplot(aes(y=counts,x=allele)) + geom_violin()


phi = 40
a = rnegbin(n=5000,mu=40,theta=40/phi)
b = rnegbin(n=5000,mu=80,theta=80/phi)



ag = data.frame(counts=c(a,b),allele=c(rep("A",length(a)),rep("B",length(b))))



m3 = glmmTMB::glmmTMB(ag$counts ~ ag$allele, dispformula = ~ag$allele,family= glmmTMB::nbinom2())
summary(m3)
m4 = glmmTMB::glmmTMB(ag$counts ~ ag$allele, dispformula = ~ag$allele,family= glmmTMB::nbinom1())
summary(m4)
#aaa$





  
m5 = glmmTMB::glmmTMB(ag$counts ~ ag$allele,family= glmmTMB::nbinom1())
summary(m5)
m6 = glmmTMB::glmmTMB(ag$counts ~ ag$allele,family= glmmTMB::nbinom2())
summary(m6)


cross_data$A$cis$ASE$ase_noise_nbin1 %>% filter(component == "disp")
cross_data$A$cis$ASE$ase_noise_nbin1 %>% filter(component != "fixed")


exp(8.21)

aa = nbin.model.resultsM=left_join(cross_data$A$cis$ASE$ase_noise_nbin1 %>% filter(term=='genoB' & component=='cond'),
                              cross_data$A$cis$ASE$ase_noise_nbin1 %>% filter(term=='genoB' & component=='disp'), by='gene', suffix=c('.cond', '.disp'))

(aa %>% mutate(p_adj_disp = p.adjust(p.value.disp,method="BH")) %>% filter(p_adj_disp < 0.05)) %>% ggplot(aes(y=estimate.disp/estimate.cond,x=estimate.cond)) + geom_point() + stat_cor()



aa



aaa = cross_data$A$cis$ASE$ase_noise %>% inner_join(aa,by="gene") 


plot(aaa$BIC.cond.x- aaa$BIC.cond.y)
abline(0,1)
#(aa %>% mutate(p_adj_disp = p.adjust(p.value.disp,method="BH")) %>% filter(p_adj_disp < 0.05)) %>% ggplot(aes(x=estimate.disp,y=estimate.cond)) + geom_point() + stat_cor()


aaa %>% filter(gene == "YFL014W")

plot(aa$estimate.disp,aa$sigma.disp)
# if likelihood ratio test
nbin.model.resultsM=nbin.model.resultsM %>% mutate(p.llrt.disp=pchisq(-2*(logLik.red.disp-logLik.disp),1,lower.tail=F))


