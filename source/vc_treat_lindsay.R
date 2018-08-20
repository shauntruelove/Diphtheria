library(grid)
library(gridExtra)
library(cowplot)
library(dplyr)

require(pracma)
alpha=0.9
Vc.Treated<-function(delay, alpha){
  load(file = "results/carriage_curve_shared.RData") # clearance_curve
  load(file = "results/carriage_curve.RData") # clearance_curve
  Treat_g1<-clearance_curve_shared %>% filter(q<=delay)
  #Treat_g1<-Treat_g1[-1,]
  Treat_g2<-clearance_curve %>% filter(type=="Case, Antibiotics")
  Treat_g2<-Treat_g2[-1,]
  Treat_g2$vcpmid<- Treat_g1$vcpmid[which(Treat_g1$q == max(Treat_g1$q))]*Treat_g2$vcpmid
  Treat_g2$vcplow<- Treat_g1$vcplow[which(Treat_g1$q == max(Treat_g1$q))]*Treat_g2$vcplow
  Treat_g2$vcphigh<- Treat_g1$vcphigh[which(Treat_g1$q == max(Treat_g1$q))]*Treat_g2$vcphigh
  Treat_g2$q<- Treat_g1$q[which(Treat_g1$q == max(Treat_g1$q))]+Treat_g2$q
  Treat_g<-rbind(Treat_g1,Treat_g2)
  Treat_g$type<-"Antibiotics"
  Case<-clearance_curve_shared %>% filter(type=="Shared") 
  #Case<-Case[-1,]
  
  (p2<- ggplot() +
      #geom_ribbon(aes(ymin=vcplow*100, ymax=vcphigh*100),  alpha=.35) +
      geom_area(data=Case, aes(x=q,y=vcpmid*100), fill="grey90", colour="grey90", linetype="dashed") +
      geom_line(data=Treat_g, aes(x=q,y=vcpmid*100), colour="Maroon", size=1) + 
      geom_area(data=Treat_g1, aes(x=q,y=vcpmid*100), colour="Maroon",fill="Maroon")+
      geom_area(data=Treat_g2, aes(x=q,y=vcpmid*100), colour="Maroon",fill="Maroon")+
      scale_x_continuous(limits=c(0, 70), expand = c(0, 0)) + 
      ylab("Percent Colonized (%)") + xlab("Duration (days)") + 
      theme(legend.position = c(0.6, 0.8),aspect.ratio = 0.6) +
      background_grid(major = "xy", minor = "none") +
      theme(text=element_text(family="Myriad Pro")))
  
  # Case<-Case %>% filter(q<=30) ## Carriage curve for untreated g(t)
  # Treat_g<-Treat_g %>% filter(q<=30) ## Carriage curve for treated g_T(T)
  ## Find the area under the curve from t=0 to t=30
  int_g_untreat<-trapz(Case$q, Case$vcpmid) 
  int_g_treat<-trapz(Treat_g$q, Treat_g$vcpmid)
  
  delta=(int_g_treat)/int_g_untreat
  
  load(file='results/R0_res_knownGT_2000iters.RData')  ## Extract data from stan
  fits <- rstan::extract(fit_knGT_fixCoV)
  
  R095 <- quantile(fits$R0_preds, probs=c(0,1))
  fits95 <- data.frame(R0_preds=fits$R0_preds, lRs_star=fits$lRs_star, Rs=fits$Rs_tmp, pk=fits$pk, 
                       tau=fits$tau, Vc=fits$Vc)
  fits95$R0 <- fits95$pk.1*fits95$Rs + (1-fits95$pk.1)*fits95$tau*fits95$Rs
  fits95 <- fits95 %>% filter(R0>=R095[1] & R0<=R095[2])
  # Calculate treatment-adjusted Vc for each iteration
  fits95$Vc.treat <- NA
  fits95$Vc <- NA
  pk.1<-fits95$pk.1
  pk.3<-fits95$pk.3
  Rs<-fits95$Rs
  tau<-fits95$tau
  for (i in 1:nrow(fits95)) {
    fits95$Vc.treat[i] <- ((-1) + pk.1[i]*Rs[i] - pk.1[i]*alpha*Rs[i] + pk.1[i]*alpha*delta*Rs[i] + tau[i]*Rs[i] - pk.1[i]*tau[i]*Rs[i])/((pk.1[i] - pk.3[i])*Rs[i]*(1 - alpha + alpha*delta - tau[i]))
    fits95$Vc[i]<-((1 - pk.1[i]*Rs[i] - Rs[i]*tau[i] + pk.1[i]*Rs[i]*tau[i])/((pk.1[i] - pk.3[i])*Rs[i]*(-1 + tau[i])))
  }
  return(fits95)
}

Vc.Treated1<-Vc.Treated(delay = 1, alpha = 0.9)
summary(Vc.Treated1$Vc.treat)
quantile(Vc.Treated1$Vc.treat, probs=c(0.025,0.975))

Vc.Treated2<-Vc.Treated(delay = 2, alpha = 0.5)
summary(Vc.Treated2$Vc.treat)
quantile(Vc.Treated2$Vc.treat, probs=c(0.025,0.975))

Vc.Treated5<-Vc.Treated(delay = 5, alpha = 0.25)
summary(Vc.Treated5$Vc.treat)
quantile(Vc.Treated5$Vc.treat, probs=c(0.025,0.975))


Vc.Treated0<-Vc.Treated(delay = 7, alpha = 1)
summary(Vc.Treated0$Vc.treat)
quantile(Vc.Treated0$Vc.treat, probs=c(0.025,0.975))
