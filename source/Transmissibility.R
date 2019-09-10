if(!require('pracma')) install.packages('pracma'); library(pracma)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('rstan')) install.packages('rstan'); library(rstan)

# Calculating Tau ---------------------------------------------------------
## Solve for: 
## X = beta_S ∫_0^30 g_S(t)dt / beta_A ∫_0^30 g_A(t)dt
## But since ∫_0^30 g_S(t)dt = ∫_0^30 g_A(t)dt, tau = 1/X
## Total number cases followed [Doull & Lara]
Sympt.cont<-758
Asympt.cont<-779

## Number of secondary [Doull & Lara]
Sympt.sec<-59
Asympt.sec<-14

sympt_second<-Sympt.sec/Sympt.cont
asympt_second<-Asympt.sec/Asympt.cont 


## Calculate the ratio of the number of symptomatic secondary infections per case to the number of symptomatic secondary infections per carrier
X= sympt_second/asympt_second

## Solve for tau = beta_A/beta_S
tau=1/(X)
tau


#Stan running
stan.data <- list(Sympt_cont= Sympt.cont,
                  Asympt_cont = Asympt.cont,
                  Sympt_sec = Sympt.sec,
                  Asympt_sec = Asympt.sec)

tau_stan <- stan("source/Transmissibility.stan", data = stan.data)
tmp <- rstan::extract(tau_stan)

#overididing old tau.df 
tau.df <- c(quantile(tmp$tau, probs=c(0.5, 0.025, 0.975)))

hist(tmp$tau, breaks=100)
mean(tau.df)
sd(tau.df)


tau.tmp <- tmp$tau
hist(log(tmp$tau), breaks=100)

logtau <- log(tmp$tau)
mean(logtau)
sd(logtau)
