load("~/Documents/diphtheria/trunk/Published Analysis/results/R0_res_knownGT_1000iters.RData")
fits <- rstan::extract(fit_knGT)
Vc<-fits$Vc
tau<-fits$tau
df<-data.frame(tau,Vc)


ggplot(df, aes(x=tau, y=Vc*100)) + 
  geom_point() +
  geom_hline(yintercept=100, color="red") + 
  ylab("Critical Vaccination Threshold (%)") +
  xlab(expression(tau))
