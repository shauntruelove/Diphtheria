# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)
if(!require('knitr')) install.packages('knitr'); library(knitr)

if(!require('fitdistrplus')) install.packages('fitdistrplus'); library(fitdistrplus)
if(!require('coda')) install.packages('coda'); library(coda)
if(!require('rjags')) install.packages('rjags'); library(rjags)

library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(knitr)

library(R2jags)
library(matrixStats)
library(loo)
library(rstanarm)
options(mc.cores = 4)



# FUNCTIONS FOR JAGS FITTING ----------------------------------------------

# Function for setting up JAGS
carriage.jags.setup.fn <- function(data, combo.data=FALSE) {
  
  #assuming observations with no upper limit have NAs in the R slot. Currently only in carr.max column
  data$carr.max[which(is.infinite(data$carr.max))] <-  NA
  
  #missing observation are just unbounded, replace with 0
  data$carr.min[which(is.na(data$carr.min))] <-  0
  
  #make a data object for JAGS
  jags.data <- list(
    CL=data$carr.min,
    CR=data$carr.max)
  
  if(combo.data) {
    jags.data <- list(
      CL=data$carr.min,
      CR=data$carr.max,
      Sympt=data$Case)
  }
  
  #Let jags know the censoring situation
  #1 means interval censored
  #2 means event occured after known time
  jags.data$CDisCensored=rep(1, length(jags.data$CL))
  jags.data$CDisCensored[which(is.na(jags.data$CR))]=2
  
  #define variable to hold the length of the time to event for bacterial clearance
  jags.data$Y_V <- rep(NA, nrow(data))
  
  #set the intitial values for time to event (i.e., Y_V)
  CDyInit =  jags.data$CL
  CDyInit[which(is.na(jags.data$CL)==T)]=0
  CDyInit[which(CDyInit==0)]=0.000000001
  
  
  #identify which observations are interval censored and which are right censored. 
  jags.data$CDic=which(jags.data$CDisCensored==1)
  jags.data$CDrc=which(jags.data$CDisCensored==2)
  
  
  if(combo.data){
    jags.data$CDic=which(jags.data$CDisCensored==1 & jags.data$Sympt==1)
    jags.data$CDrc=which(jags.data$CDisCensored==2 & jags.data$Sympt==1)
    jags.data$CDic2=which(jags.data$CDisCensored==1 & jags.data$Sympt==0)
    jags.data$CDrc2=which(jags.data$CDisCensored==2 & jags.data$Sympt==0)
  }
  
  # add to right bound to fix problem with having left and right bounds the same
  jags.data$CR <- jags.data$CR + (23/24)
  
  #set the parameters we want to track
  #  - v_v     : weibull shape
  #  - scale_v : weibull scale
  parameters <- c("v_v","scale_v", "loglik")  
  
  set.seed(12345) #if this is not included, multiple initializations maybe needed.
  #initialization function for jags
  return(list(jags.data=jags.data, parameters=parameters))
}


# LOAD & SETUP DATA ---------------------------------------------------------------

# load the raw colonization data
dc.data <- read.csv(file = "./data/duration_colonization.csv", header=TRUE)
#dc.data.treat <- read.csv(file = "./data/duration_colonization_treated.csv", header=TRUE)

# Convert to long
dc.long <- data.frame(Case=rep(dc.data$Case, times=dc.data$Num.convert),
                      Treatment=rep(dc.data$Treatment, times=dc.data$Num.convert),
                      Antitoxin=rep(dc.data$Antitoxin, times=dc.data$Num.convert),
                      carr.min=rep(dc.data$carr.min, times=dc.data$Num.convert), 
                      carr.max=rep(dc.data$carr.max, times=dc.data$Num.convert),
                      study=rep(dc.data$Study, times=dc.data$Num.convert))
dc.long$Case <- ifelse(dc.long$Case=="Case", 1, 0)
dc.long <- dc.long %>% filter(Treatment=="none")
dc.long.sympt <- dc.long %>% filter(Case==1 & Treatment=="none")
dc.long.asympt <- dc.long %>% filter(Case==0 & Treatment=="none")

# Setup data
jags.setup.shared <- carriage.jags.setup.fn(dc.long)
jags.setup.asympt <- carriage.jags.setup.fn(dc.long.asympt)
jags.setup.sympt <- carriage.jags.setup.fn(dc.long.sympt)

jags.setup.comb <- carriage.jags.setup.fn(dc.long, combo.data = TRUE)


# RUN JAGS MODEL FOR EACH -------------------------------------------------

# Shared fit
shared.fit <- jags(data = jags.setup.shared$jags.data,
                    parameters.to.save = jags.setup.shared$parameters,
                    n.chains = 3,
                    n.iter = 10000,
                    n.burnin = 1000,
                    n.thin = 3,
                    model.file = './source/DistributionFitCarriage.jags')

# Symptomatic fit
sympt.fit <- jags(data = jags.setup.sympt$jags.data,
                  parameters.to.save = jags.setup.sympt$parameters,
                  n.chains = 3,
                  n.iter = 10000,
                  n.burnin = 1000,
                  n.thin = 3,
                  model.file = './source/DistributionFitCarriage.jags')

# Asymptomatic fit
asympt.fit <- jags(data = jags.setup.asympt$jags.data,
                  parameters.to.save = jags.setup.asympt$parameters,
                  n.chains = 3,
                  n.iter = 10000,
                  n.burnin = 1000,
                  n.thin = 3,
                  model.file = './source/DistributionFitCarriage.jags')

# Combo fit
comb.fit <- jags(data = jags.setup.comb$jags.data,
                   parameters.to.save = jags.setup.comb$parameters,
                   n.chains = 3,
                   n.iter = 10000,
                   n.burnin = 1000,
                   n.thin = 3,
                   model.file = './source/DistributionFitCarriage_combined.jags')



# CHECK WAIC between datasets ---------------------------------------------

shared.paramlist <- shared.fit$BUGSoutput$sims.list
shared.loglik <- shared.paramlist$loglik

sympt.paramlist <- sympt.fit$BUGSoutput$sims.list
sympt.loglik <- sympt.paramlist$loglik

asympt.paramlist <- asympt.fit$BUGSoutput$sims.list
asympt.loglik <- asympt.paramlist$loglik

comb.paramlist <- comb.fit$BUGSoutput$sims.list
comb.loglik <- comb.paramlist$loglik



# Combine the loglikelihood matrices for asympt and sympt fits.
comb.loglik2 <- cbind(sympt.loglik, asympt.loglik)

# WAIC from loo package
shared.waic <- waic(shared.loglik)
comb.waic <- waic(comb.loglik)
comb.waic2 <- waic(comb.loglik2)

# LOO from loo package
shared.loo <- loo(shared.loglik)
comb.loo <- loo(comb.loglik)

# Compare results from fits
compare(shared.waic, comb.waic)
compare(comb.waic, comb.waic2)
compare(shared.loo, comb.loo)
compare(shared.waic, comb.waic2)


# Shared distribution does a better/equivalent job of fitting 
#   asymptomatic and symptomatic colonization duration.






