# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)



# Proportion who are asymptomatic by vaccination status -----------------------

Asympt <- read.csv("data/asymptomatic_by_vaccination.csv", header = TRUE, stringsAsFactors = FALSE)
Asympt <- data.frame(Asympt)


# For only non-vaccinated  ------------------------------------------------
asym <- sum(Asympt$carr_non)
total <- sum(Asympt$carr_non,Asympt$case_non)

prop_Asympt_non<- 1- asym/total

prop_Asympt_non_LCI <- prop_Asympt_non - 1.96*sqrt(prop_Asympt_non*(1-prop_Asympt_non)/total)
prop_Asympt_non_UCI <- prop_Asympt_non + 1.96*sqrt(prop_Asympt_non*(1-prop_Asympt_non)/total)



# Combine into table ------------------------------------------------------

Proportion_Asymptomatic <- data.frame(prop_Asympt_non, prop_Asympt_non_LCI, prop_Asympt_non_UCI, "P0")
colnames(Proportion_Asymptomatic) <- c("Proportion", "Lower CI", "Upper CI", "Type")

# Proportion Asymptomatic among those not vaccinated
print(Proportion_Asymptomatic)





# Bring in Stan results to estimate p2 and p3 -----------------------------
#  --- From the Stan results we can extract distributions for p{1-2} and p3

load(file='results/R0_res_knownGT_2000iters.RData') # Loads ret_sm, fit_knGT, full_data
chains <- rstan::extract(fit_knGT)

p1 <- chains$pk[,2]
p3 <- chains$pk[,3]

summary(p1)
quantile(p1,probs=c(0.025,0.975))

summary(p3)
quantile(p3,probs=c(0.025,0.975))

##





















# # Relative Risk ------------------------------------------------------------
# head(Asympt)
# Asympt<-data.frame(Asympt)
# Asympt2<-Asympt[-c(1,4),]
# 
# TotAsympt<-matrix(NA,2,2)
# TotAsympt[1,1]<-sum(Asympt2$case_full)
# TotAsympt[2,1]<-sum(Asympt2$case_partial,Asympt2$case_non)
# TotAsympt[1,2]<-sum(Asympt2$carr_full)
# TotAsympt[2,2]<-sum(Asympt2$carr_partial,Asympt2$carr_non)
# 
# TotAsympt<-data.frame(TotAsympt)
# TotAsympt<-as.numeric(TotAsympt)
# colnames(TotAsympt)<-c( "Sympt", "Asympt")
# 
# require(fmsb)
# 
# X<-TotAsympt[1,2]
# Y<-TotAsympt[2,2]
# M1<-sum(TotAsympt[1,])
# M2<-sum(TotAsympt[2,])
# riskratio(X, Y, M1, M2, conf.level=0.95)
# 
