


calc_delta <- function(delay=1){
  
  require(pracma)
  
  # load(file = "results/carriage_curve_shared.RData") # clearance_curve_shared
  # load(file = "results/carriage_curve_treated.RData") # clearance_curve_treated
  clearance_curve_shared <- read.csv(file='results/clearance_curve_shared.csv', header=TRUE)
  clearance_curve_treated <- read.csv(file='results/clearance_curve_treated.csv', header=TRUE)
  
  Treat_g1 <- clearance_curve_shared %>% filter(q<=delay)
  Treat_g2 <- clearance_curve_treated
  Treat_g2$vcpmid<- Treat_g1$vcpmid[which(Treat_g1$q == max(Treat_g1$q))]*Treat_g2$vcpmid
  
  Treat_g2$vcplow <- Treat_g1$vcplow[which(Treat_g1$q == max(Treat_g1$q))]*Treat_g2$vcplow
  Treat_g2$vcphigh <- Treat_g1$vcphigh[which(Treat_g1$q == max(Treat_g1$q))]*Treat_g2$vcphigh
  Treat_g2$q <- Treat_g1$q[which(Treat_g1$q == max(Treat_g1$q))]+Treat_g2$q
  Treat_g <- rbind(Treat_g1,Treat_g2)
  Treat_g$type<-"Antibiotics"
  Case <- clearance_curve_shared %>% filter(type=="Shared") 
  #Case <- Case[-1,]
  
  int_g_untreat <- trapz(Case$q, Case$vcpmid) 
  int_g_treat <- trapz(Treat_g$q, Treat_g$vcpmid)
  
  delta <- (int_g_treat)/int_g_untreat
  
  return(delta)
}


vc_treat <- function(delay, alpha, restrict_quantiles=c(0,1)){
  
  delta <- calc_delta(delay=delay)
  
  load(file='results/R0_res_knownGT_2000iters.RData')  ## Extract data from stan
  fits <- rstan::extract(fit_knGT_fixCoV)

  R0_quantiles <- quantile(fits$R0_preds, probs=restrict_quantiles)
  fits_data <- data.frame(R0_preds=fits$R0_preds, lRs_star=fits$lRs_star, Rs=fits$Rs_tmp, p0=fits$pk[,1], p3=fits$pk[,3], 
                          tau=fits$tau, Vc=fits$Vc)
  fits_data$R0 <- fits_data$p0*fits_data$Rs + (1-fits_data$p0)*fits_data$tau*fits_data$Rs
  fits_data <- fits_data %>% filter(R0>=R0_quantiles[1] & R0<=R0_quantiles[2])
  
  # Calculate treatment-adjusted Vc for each iteration
  fits_data$Vc.treat <- NA
  fits_data$Vc <- NA
  p0 <- fits_data$p0
  p3 <- fits_data$p3
  Rs <- fits_data$Rs
  tau <- fits_data$tau
  for (i in 1:nrow(fits_data)) {
    #fits_data$Vc.treat[i] <- calc_Vc_treatment(p0=p0[i], p3=p3[i], Rs=Rs[i], Ra=NULL, tau=tau[i], delta=delta, alpha=alpha)
    fits_data$Vc.treat[i] <- ((-1) + p0[i]*Rs[i] - p0[i]*alpha*Rs[i] + p0[i]*alpha*delta*Rs[i] + tau[i]*Rs[i] - p0[i]*tau[i]*Rs[i])/((p0[i] - p3[i])*Rs[i]*(1 - alpha + alpha*delta - tau[i]))
    fits_data$Vc[i]<-((1 - p0[i]*Rs[i] - Rs[i]*tau[i] + p0[i]*Rs[i]*tau[i])/((p0[i] - p3[i])*Rs[i]*(-1 + tau[i])))
  }
  return(fits_data)
}


# Function to calc Vc with treatment
calc_Vc_treatment <- function(p0, p3, Rs, Ra=NULL, tau, delta, alpha){
  if (is.null(Ra))  Ra <- tau*Rs
  # return((1/Rs - (p0*(1-alpha) + p0*alpha*delta + tau - tau*p0)) /
  #            (-p0*(1-alpha) - p0*alpha*delta - tau + tau*p0 + p3*(1-alpha) + p3*alpha*delta + tau - p3*tau))
  return((-1 - p0*Rs - p0*alpha*Rs + p0*alpha*delta*Rs - Rs*tau - p0*Rs*tau) / ((p0 - p3)*Rs*(1 - alpha + alpha*tau - tau)))
}


# Function to calc Vc
calc_Vc <- function(p0, p3, Rs, Ra=NULL, tau){
  if (is.null(Ra))  Ra <- tau*Rs
  return((1 - p0*Rs - Ra + p0*Ra) / (p3*Rs - p3*Ra - p0*Rs + p0*Ra))
}


# Function to estimate the overall adjusted Rs in the population (adjusted for antibiotics)
Rs_treat <- function(delta, alpha, Rs) {
  delta*alpha*Rs + (1-alpha)*Rs
}

# Proportion to treat to drop R<1
calc_alpha <- function(p3, Rs, tau, delta){
  return((1 - p3*Rs - Rs*tau + p3*Rs*tau)/((-1 + delta)*p3*Rs))
}

# Proportion to treat to drop R<1 -- with Ra treatment
calc_alpha_RsRa <- function(p3, Rs, tau, delta){
  return((-1 + Rs*(p3 + tau - p3*tau))/((-1 + delta)*Rs*(p3*(-1 + tau) - tau)))
}


# Set up Functions --------------------------------------------------------

Vc_funct <- function(Rs, p0, p3, tau, delta, alpha){
  ((1/Rs) - alpha*p0 - delta*p0 + alpha*delta*p0 - tau + p0*tau) / 
    ((p0-p3)*((-alpha) - delta + alpha*delta + tau))
}

R_calc <- function(V3,p0,p3,Rs,tau,alpha,delta){
  Rs_star <- alpha*delta*Rs + (1-alpha)*Rs
  Ra <- tau*Rs
  (1-V3)*(p0*Rs_star + (1-p0)*Ra) + V3*(p3*Rs_star + (1-p3)*Ra)
}

# Proportion to treat to drop R<1
calc_alpha <- function(p3, Rs, tau, delta){
  return((1 - p3*Rs - Rs*tau + p3*Rs*tau)/((-1 + delta)*p3*Rs))
}

get_R0_IQR <- function(delay, alpha, restrict_quantiles=c(0,1)){
  R0.data <- vc_treat(delay, alpha, restrict_quantiles)$R0
  return(quantile(R0.data, probs=c(0.25,0.75)))
}





# 
# 
# vc_treat_R0 <- function(delay, alpha, restrict_quantiles=c(0,1)){
#   
#   load(file = "results/carriage_curve_shared.RData") # clearance_curve
#   load(file = "results/carriage_curve.RData") # clearance_curve
#   
#   Treat_g1 <- clearance_curve_shared %>% filter(q<=delay)
#   Treat_g1 <- Treat_g1[-1,]
#   Treat_g2 <- clearance_curve %>% filter(type=="Case, Antibiotics")
#   Treat_g2 <- Treat_g2[-1,]
#   Treat_g2$vcpmid<- Treat_g1$vcpmid[which(Treat_g1$q == max(Treat_g1$q))]*Treat_g2$vcpmid
#   
#   Treat_g2$vcplow <- Treat_g1$vcplow[which(Treat_g1$q == max(Treat_g1$q))]*Treat_g2$vcplow
#   Treat_g2$vcphigh <- Treat_g1$vcphigh[which(Treat_g1$q == max(Treat_g1$q))]*Treat_g2$vcphigh
#   Treat_g2$q <- Treat_g1$q[which(Treat_g1$q == max(Treat_g1$q))]+Treat_g2$q
#   Treat_g <- rbind(Treat_g1,Treat_g2)
#   Treat_g$type<-"Antibiotics"
#   Case <- clearance_curve_shared %>% filter(type=="Shared") 
#   Case <- Case[-1,]
#   
#   (p2<- ggplot() +
#       #geom_ribbon(aes(ymin=vcplow*100, ymax=vcphigh*100),  alpha=.35) +
#       geom_area(data=Case, aes(x=q,y=vcpmid*100), fill="grey90", colour="grey90", linetype="dashed") +
#       geom_line(data=Treat_g, aes(x=q,y=vcpmid*100), colour="Maroon", size=1) + 
#       geom_area(data=Treat_g1, aes(x=q,y=vcpmid*100), colour="Maroon",fill="Maroon")+
#       geom_area(data=Treat_g2, aes(x=q,y=vcpmid*100), colour="Maroon",fill="Maroon")+
#       scale_x_continuous(limits=c(0, 70), expand = c(0, 0)) + 
#       ylab("Percent Colonized (%)") + xlab("Duration (days)") + 
#       theme(legend.position = c(0.6, 0.8),aspect.ratio = 0.6) +
#       background_grid(major = "xy", minor = "none") +
#       theme(text=element_text(family="Myriad Pro")))
#   
#   # Case <- Case %>% filter(q<=30) ## Carriage curve for untreated g(t)
#   # Treat_g <- Treat_g %>% filter(q<=30) ## Carriage curve for treated g_T(T)
#   ## Find the area under the curve from t=0 to t=30
#   int_g_untreat <- trapz(Case$q, Case$vcpmid) 
#   int_g_treat <- trapz(Treat_g$q, Treat_g$vcpmid)
#   
#   delta=(int_g_treat)/int_g_untreat
#   
#   load(file='results/R0_res_knownGT_2000iters.RData')  ## Extract data from stan
#   fits <- rstan::extract(fit_knGT_fixCoV)
#   
#   R0_quantiles <- quantile(fits$R0_preds, probs=restrict_quantiles)
#   fits_data <- data.frame(R0_preds=fits$R0_preds, lRs_star=fits$lRs_star, Rs=fits$Rs_tmp, pk=fits$pk, 
#                           tau=fits$tau, Vc=fits$Vc)
#   fits_data$R0 <- fits_data$p0*fits_data$Rs + (1-fits_data$p0)*fits_data$tau*fits_data$Rs
#   fits_data <- fits_data %>% filter(R0>=R0_quantiles[1] & R0<=R0_quantiles[2])
#   # Calculate treatment-adjusted Vc for each iteration
#   fits_data$Vc.treat <- NA
#   fits_data$Vc <- NA
#   p0 <- fits_data$p0
#   p3 <- fits_data$p3
#   Rs <- fits_data$Rs
#   tau <- fits_data$tau
#   for (i in 1:nrow(fits_data)) {
#     fits_data$Vc.treat[i] <- ((-1) + p0[i]*Rs[i] - p0[i]*alpha*Rs[i] + p0[i]*alpha*delta*Rs[i] + tau[i]*Rs[i] - p0[i]*tau[i]*Rs[i])/((p0[i] - p3[i])*Rs[i]*(1 - alpha + alpha*delta - tau[i]))
#     fits_data$Vc[i]<-((1 - p0[i]*Rs[i] - Rs[i]*tau[i] + p0[i]*Rs[i]*tau[i])/((p0[i] - p3[i])*Rs[i]*(-1 + tau[i])))
#   }
#   return(fits_data)
# }
