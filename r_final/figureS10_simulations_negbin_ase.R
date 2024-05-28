

p1 = neg_bin_sims$sim1 %>% ggplot(aes(x=-theta2,y=-delta.disp, group=factor(fc), color=factor(fc))) + geom_point() + theme_bw() + ylab(expression("estimated "*Delta * ln(noise))) +
    xlab(expression("simulated "*Delta * ln(noise))) + ggtitle(paste("Baseline expression:", round(exp(unique(neg_bin_sims$sim1$mu)),digits = 2),"counts per cell")) +
  viridis::scale_color_viridis(name=expression("simulated "*Delta * ln(expression)),discrete = T) + 
  geom_abline(intercept = 0,slope=1) + geom_line() #+ scale_colo
p2 = neg_bin_sims$sim1 %>% ggplot(aes(x=fc,y=delta.mean, group=factor(-theta2), color=factor(-theta2))) + geom_point() + theme_bw() + ylab(expression("estimated "*Delta * ln(expression))) +
  xlab(expression("simulated "*Delta * ln(expression))) + ggtitle(paste("Baseline expression:", round(exp(unique(neg_bin_sims$sim1$mu)),digits = 2),"counts per cell")) +
  viridis::scale_color_viridis(name=expression("simulated "*Delta * ln(noise)),discrete = T) + 
  geom_abline(intercept = 0,slope=1) + geom_line() #+ scale_colo
p3 = neg_bin_sims$sim2 %>% ggplot(aes(x=-theta2,y=-delta.disp, group=factor(fc), color=factor(fc))) + geom_point() + theme_bw() + ylab(expression("estimated "*Delta * ln(noise))) +
  xlab(expression("simulated "*Delta * ln(noise))) + ggtitle(paste("Baseline expression:", round(exp(unique(neg_bin_sims$sim2$mu)),digits = 2),"counts per cell")) + 
  viridis::scale_color_viridis(name=expression("simulated "*Delta * ln(expression)),discrete = T) +
   geom_abline(intercept = 0,slope=1) + geom_line()
#p4 = neg_bin_sims$sim2 %>% ggplot(aes(x=fc,y=delta.mean, group=factor(-theta2), color=factor(-theta2))) + geom_point() + theme_bw() + ylab(expression("estimated "*Delta * ln(noise))) +
#  xlab(expression("simulated "*Delta * ln(noise))) + ggtitle(paste("Baseline expression:", round(exp(unique(neg_bin_sims$sim2$mu)),digits = 2),"counts per cell")) +
#  viridis::scale_color_viridis(name="ln(noise)",discrete = T) + 
#  geom_abline(intercept = 0,slope=1) + geom_line()

p4 = neg_bin_sims$sim2 %>% ggplot(aes(x=fc,y=delta.mean, group=factor(-theta2), color=factor(-theta2))) + geom_point() + theme_bw() + ylab(expression("estimated "*Delta * ln(expression))) +
  xlab(expression("simulated "*Delta * ln(expression))) + ggtitle(paste("Baseline expression:", round(exp(unique(neg_bin_sims$sim2$mu)),digits = 2),"counts per cell")) +
  viridis::scale_color_viridis(name=expression("simulated "*Delta * ln(noise)),discrete = T) + 
  geom_abline(intercept = 0,slope=1) + geom_line() #+ scale_colo


#
#neg_bin_sims$sim2$mu

plot_grid(p1,p2,p3,p4,nrow=2,labels = "AUTO",label_size = 16)
ggsave("fig_final/s10.png",width=16,height=12,dpi=300,device = png)
ggsave("fig_final//svg/s10.svg",width=16,height=12,device = png)
