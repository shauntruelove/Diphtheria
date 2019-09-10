# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)
if(!require('fitdistrplus')) install.packages('fitdistrplus'); library(fitdistrplus)
if(!require('coda')) install.packages('coda'); library(coda)
if(!require('rjags')) install.packages('rjags'); library(rjags)


# Source files & data ------------------------------------------------------------

source('source/Vc_and_R0_functions.R')

load(file='results/R0_res_knownGT_2000iters.RData')  ## Extract data from stan
fits <- rstan::extract(fit_knGT_fixCoV)

rstan::summary(fit_knGT_fixCoV)

# Estimate relative infectiousness of treated individuals (delta) -------------------------------------

delta1day <- calc_delta(delay=1)
delta2day <- calc_delta(delay=2)
delta5day <- calc_delta(delay=5)


# Clean Data and Set Up as 95% Dataframe ----------------------------------

# Get the stan results and restrict to the 95% CI of R0

R095 <- quantile(fits$R0_preds, probs=c(0.025,0.975))
fits95 <- data.frame(R0_preds=fits$R0_preds, lRs_star=fits$lRs_star, Rs=fits$Rs_tmp, 
                     p0=fits$pk[,1], p2=fits$pk[,2], p3=fits$pk[,3],
                     tau=fits$tau, Vc=fits$Vc)
fits95$R0 <- fits95$p0*fits95$Rs + (1-fits95$p0)*fits95$tau*fits95$Rs # Double Check R0
fits95 <- fits95 %>% filter(R0>=R095[1] & R0<=R095[2])


#vc_treat <- function(delay, alpha, restrict_quantiles=c(0,1)){
  



# Calculate various metrics -----------------------------------------------


# Calculate treatment-adjusted Vc for each iteration

# 3 Scenarios
fits95$Vc_scen1 <- fits95$Vc_scen2 <- fits95$Vc_scen3 <- NA
alpha <- c(.9, .5, .25)
delay <- c(1,2,5)
for (i in 1:nrow(fits95)) {
    Rs_star_scen1 <- Rs_treat(delta=calc_delta(delay=delay[1]), alpha=alpha[1], Rs=fits95$Rs[i])
    Rs_star_scen2 <- Rs_treat(delta=calc_delta(delay=delay[2]), alpha=alpha[2], Rs=fits95$Rs[i])
    Rs_star_scen3 <- Rs_treat(delta=calc_delta(delay=delay[3]), alpha=alpha[3], Rs=fits95$Rs[i])
    
    fits95$Vc_scen1[i] <- calc_Vc(p0=fits95$p0[i], p3=fits95$p3[i], Rs=Rs_star_scen1, Ra=(fits95$tau[i]*fits95$Rs[i]), tau=fits95$tau[i])
    fits95$Vc_scen2[i] <- calc_Vc(p0=fits95$p0[i], p3=fits95$p3[i], Rs=Rs_star_scen2, Ra=(fits95$tau[i]*fits95$Rs[i]), tau=fits95$tau[i])
    fits95$Vc_scen3[i] <- calc_Vc(p0=fits95$p0[i], p3=fits95$p3[i], Rs=Rs_star_scen3, Ra=(fits95$tau[i]*fits95$Rs[i]), tau=fits95$tau[i])
}
    
for (i in 1:nrow(fits95)) {
    fits95$PropIsolate[i] <- 1-calc_alpha(p3=fits95$p3[i], Rs=fits95$Rs[i], tau=fits95$tau[i], delta=0)
}

for (i in 1:nrow(fits95)) {
  fits95$PropTreat_Mass[i] <- calc_alpha_RsRa(p3=fits95$p3[i], Rs=fits95$Rs[i], tau=fits95$tau[i], delta=calc_delta(delay=10-1.9))
  #fits95$PropTreat_Mass[i] <- calc_alpha_RsRa(p3=fits95$p3[i], Rs=fits95$Rs[i], tau=fits95$tau[i], delta=calc_delta(treat.delay=10, time.to.membr=0))
}




View(fits95)


hist(fits95$Vc, breaks=100)
hist(fits95$Vc_scen1, breaks=100)
hist(fits95$Vc_scen2, breaks=100)
hist(fits95$Vc_scen3, breaks=100)

summary(fits95$Vc)
# summary(fits95$Vc_adj)
quantile(fits95$Vc, probs=c(.025,0.975))


summary(fits95$PropIsolate)
quantile(fits95$PropIsolate, probs=c(0.025,0.975))


# Proportion of Vc>1 ------------------------------------------------------

mean(fits95$Vc<1)
mean(fits95$Vc_scen1<1)
mean(fits95$Vc_scen2<1)
mean(fits95$Vc_scen3<1)

median(fits95$Vc)
median(fits95$Vc_scen1)
median(fits95$Vc_scen2)
median(fits95$Vc_scen3)

mean(fits95$Vc)
mean(fits95$Vc_scen1)
mean(fits95$Vc_scen2)
mean(fits95$Vc_scen3)

quantile(fits95$Vc_scen1, probs=c(0.025, .975))
quantile(fits95$Vc_scen2, probs=c(0.025, .975))
quantile(fits95$Vc_scen3, probs=c(0.025, .975))



# Mass Antibiotics --------------------------------------------------------

summary(fits95$PropTreat_Mass)
quantile(fits95$PropTreat_Mass, probs=c(0.025,0.975))




# Proportion Need to Treat to have Vc<100% and R<1 ------------------------

fits95$alpha_delta1day <- NA
fits95$alpha_delta2day <- NA
fits95$alpha_delta5day <- NA

for (i in 1:nrow(fits95)) {
    fits95$alpha_delta1day[i] <- calc_alpha(fits95$p3[i], fits95$Rs[i], fits95$tau[i], calc_delta(treat.delay=1, time.to.membr=(1.9)))
    fits95$alpha_delta2day[i] <- calc_alpha(fits95$p3[i], fits95$Rs[i], fits95$tau[i], calc_delta(treat.delay=2, time.to.membr=(1.9)))
    fits95$alpha_delta5day[i] <- calc_alpha(fits95$p3[i], fits95$Rs[i], fits95$tau[i], calc_delta(treat.delay=5, time.to.membr=(1.9)))
}


plot(fits95$R0, fits95$alpha_delta1day, ylim=c(-1,2))

mean(fits95$alpha_delta1day<1)
mean(fits95$alpha_delta2day<1)
mean(fits95$alpha_delta5day<1)

mean(fits95$alpha_delta2day)
median(fits95$alpha_delta2day)
quantile(fits95$alpha_delta2day, probs = c(0.025,.975))

mean(fits95$alpha_delta2day)
median(fits95$alpha_delta2day)
quantile(fits95$alpha_delta2day, probs = c(0.025,.975))















# Run 3 treatment scenarios (25, 50, 75% antibiotic treatment) ------------


# Treatment rate (with antibiotics)
alpha=c(.25, 0.5, 0.75, .90)

# Get the stan results and restrict to the 95% CI of R0
R095 <- quantile(fits$R0_preds, probs=c(0.025,0.975))
Vc_treat_1day <- data.frame(R0_preds=fits$R0_preds, lRs_star=fits$lRs_star, Rs=fits$Rs_tmp, pk=fits$pk, tau=fits$tau, Vc=fits$Vc)
Vc_treat_1day$R0 <- Vc_treat_1day$pk.1*Vc_treat_1day$Rs + (1-Vc_treat_1day$pk.1)*Vc_treat_1day$tau*Vc_treat_1day$Rs
Vc_treat_1day <- Vc_treat_1day %>% filter(R0>=R095[1] & R0<=R095[2])

# Calculate treatment-adjusted Vc for each iteration
Vc_treat_1day$Vc_25pct <- Vc_treat_1day$Vc_50pct <- Vc_treat_1day$Vc_75pct <- NA
Vc_treat_1day$Rs_25pct <- Vc_treat_1day$Rs_50pct <- Vc_treat_1day$Rs_75pct <- NA


for (i in 1:nrow(Vc_treat_1day)) {
    # Vc_treat_1day$Rs_25pct[i] <- Rs_treat(delta=prop.redux.2day, alpha=.25, Rs=Vc_treat_1day$Rs[i])
    # Vc_treat_1day$Rs_50pct[i] <- Rs_treat(delta=prop.redux.2day, alpha=.50, Rs=Vc_treat_1day$Rs[i])
    # Vc_treat_1day$Rs_75pct[i] <- Rs_treat(delta=prop.redux.2day, alpha=.75, Rs=Vc_treat_1day$Rs[i])
    calc_Vc(p0=Vc_treat_1day$pk.1[i], p3=Vc_treat_1day$pk.3[i], 
                      Rs=Vc_treat_1day$Rs[i], Ra=(Vc_treat_1day$tau[i]*Vc_treat_1day$Rs[i]), 
                      tau=Vc_treat_1day$tau[i])
    
    calc_Vc_treatment(p0=Vc_treat_1day$pk.1[i], p3=Vc_treat_1day$pk.3[i], 
                      Rs=Vc_treat_1day$Rs[i], Ra=(Vc_treat_1day$tau[i]*Vc_treat_1day$Rs[i]), 
                      tau=Vc_treat_1day$tau[i], delta=seq(0,1,.1), alpha=1)
        
    Vc_treat_1day$Vc_25pct[i] <- calc_Vc(p0=Vc_treat_1day$pk.1[i], p3=Vc_treat_1day$pk.3[i], Rs=Vc_treat_1day$Rs_25pct[i], Ra=(Vc_treat_1day$tau[i]*Vc_treat_1day$Rs[i]), tau=Vc_treat_1day$tau[i])
    Vc_treat_1day$Vc_50pct[i] <- calc_Vc(p0=Vc_treat_1day$pk.1[i], p3=Vc_treat_1day$pk.3[i], Rs=Vc_treat_1day$Rs_50pct[i], Ra=(Vc_treat_1day$tau[i]*Vc_treat_1day$Rs[i]), tau=Vc_treat_1day$tau[i])
    Vc_treat_1day$Vc_75pct[i] <- calc_Vc(p0=Vc_treat_1day$pk.1[i], p3=Vc_treat_1day$pk.3[i], Rs=Vc_treat_1day$Rs_75pct[i], Ra=(Vc_treat_1day$tau[i]*Vc_treat_1day$Rs[i]), tau=Vc_treat_1day$tau[i])
}


p0=Vc_treat_1day$pk.1[i]; p3=Vc_treat_1day$pk.3[i]; 
    Rs=Vc_treat_1day$Rs[i]; Ra=(Vc_treat_1day$tau[i]*Vc_treat_1day$Rs[i]); 
    tau=Vc_treat_1day$tau[i]; delta=seq(0,1,.1); alpha=1
((1/Rs) - (p0*delta + (1-p0)*tau)) / (-p0*delta - tau*(1-p0) + p3*delta + (1-p3)*tau)

delta <- seq(0,1,.0001)
y=calc_Vc(p0=Vc_treat_1day$pk.1[i], p3=Vc_treat_1day$pk.3[i], Rs=Rs_treat(delta=delta, alpha=alpha, Rs=Vc_treat_1day$Rs[i]),
          Ra=(Vc_treat_1day$tau[i]*Vc_treat_1day$Rs[i]), tau=Vc_treat_1day$tau[i])
plot(x=delta, y)


Vc_treat_3day <- Vc_treat_1day

# Calculate treatment-adjusted Vc for each iteration
Vc_treat_3day$Vc_25pct <- Vc_treat_3day$Vc_50pct <- Vc_treat_3day$Vc_75pct <- NA
Vc_treat_3day$Rs_25pct <- Vc_treat_3day$Rs_50pct <- Vc_treat_3day$Rs_75pct <- NA


for (i in 1:nrow(Vc_treat_3day)) {
    Vc_treat_1day$Rs_25pct[i] <- Rs_treat(delta=delta1day, alpha=.25, Rs=Vc_treat_3day$Rs[i])
    Vc_treat_1day$Rs_50pct[i] <- Rs_treat(delta=delta1day, alpha=.50, Rs=Vc_treat_3day$Rs[i])
    Vc_treat_1day$Rs_75pct[i] <- Rs_treat(delta=delta1day, alpha=.75, Rs=Vc_treat_3day$Rs[i])
    
    Vc_treat_1day$Vc_25pct[i] <- calc_Vc(p0=Vc_treat_1day$pk.1[i], p3=Vc_treat_1day$pk.3[i], Rs=Vc_treat_1day$Rs_25pct[i], Ra=(Vc_treat_1day$tau[i]*Vc_treat_1day$Rs[i]), tau=Vc_treat_1day$tau[i])
    Vc_treat_1day$Vc_50pct[i] <- calc_Vc(p0=Vc_treat_1day$pk.1[i], p3=Vc_treat_1day$pk.3[i], Rs=Vc_treat_1day$Rs_50pct[i], Ra=(Vc_treat_1day$tau[i]*Vc_treat_1day$Rs[i]), tau=Vc_treat_1day$tau[i])
    Vc_treat_1day$Vc_75pct[i] <- calc_Vc(p0=Vc_treat_1day$pk.1[i], p3=Vc_treat_1day$pk.3[i], Rs=Vc_treat_1day$Rs_75pct[i], Ra=(Vc_treat_1day$tau[i]*Vc_treat_1day$Rs[i]), tau=Vc_treat_1day$tau[i])
}



save(Vc_treat_3day, Vc_treat_2day, file="./trunk/results/Vc_treat_2&3days.RData")




    

