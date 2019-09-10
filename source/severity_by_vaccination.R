library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(knitr)

# Install packages
if(!require('googlesheets')) install.packages('googlesheets'); library(googlesheets)
if(!require('fitdistrplus')) install.packages('fitdistrplus'); library(fitdistrplus)
if(!require('coda')) install.packages('coda'); library(coda)
if(!require('rjags')) install.packages('rjags'); library(rjags)



# Load and Clean Data ------------------------------------------

Vax_Sever <- read.csv(file="data/severity_by_vaccination.csv", header = TRUE, stringsAsFactors = FALSE)

Vax_Sever$Severe[is.na(Vax_Sever$Severe)] <- 0
Vax_Sever$Moderate[is.na(Vax_Sever$Moderate)] <- 0
Vax_Sever$Not.Severe..mild.[is.na(Vax_Sever$Not.Severe..mild.)] <- 0
Vax_Sever$N_group <- Vax_Sever$Severe + Vax_Sever$Moderate + Vax_Sever$Not.Severe..mild.

Vax_Sever$Population <- factor(paste(Vax_Sever$Location, Vax_Sever$Study, sep=" "))
Vax_Sever$Immunization.Status <- factor(Vax_Sever$Immunization.Status, levels=c("None","Imperfect","Full"))
Vax_Sever$N_severe <- Vax_Sever$Severe
Vax_Sever$N_notsevere <- Vax_Sever$N_group - Vax_Sever$N_severe


# Convert to long data
Vax_Sever <- Vax_Sever %>% mutate(row = 1:nrow(Vax_Sever))
Vax_Sever_long <- Vax_Sever %>% group_by(row) %>%
  do({ left_join(data_frame(row=.$row, severe = c(rep(1,.$Severe), rep(0,.$N_notsevere))),.,by='row') })




# Run GLM to get Proportions ----------------------------------------------

lm1 <- glm(severe ~ Immunization.Status, data=Vax_Sever_long, family=poisson(link='log'))
summary(lm1)

exp(coef(lm1))


# Run GLMER to get Proportions - Mixed Effects ----------------------------------------------
library(lme4)

lm2 <- glmer(severe ~ Immunization.Status + (1|Population), data=Vax_Sever_long, family=poisson(link='log'))
summary(lm2)

exp(fixef(lm2))
1-exp(fixef(lm2)) # VEs

# Relative Risks of Severity
exp(cbind("Risk ratio" = fixef(lm2), confint(lm2, method="Wald", level = 0.95)[-1,]))
1 - exp(cbind("Risk ratio" = fixef(lm2), confint(lm2, method="Wald", level = 0.95)[-1,]))

# Convert to Probabilities
library(effects)
(ef.1 <- data.frame(effect(c("Immunization.Status"), lm2)))





# 
# 
# VaxSever<-Vax_Sever[-c(2:3, 5:10,15)]
# VaxSever<- VaxSever %>% filter(Location!=c("Total"))
# VaxSever<- VaxSever %>% filter(Location!=c("Weighted"))
# 
# VaxSeverTotal<-gs_read(ss=diph.sheet, ws="Vaccination_Illness-Severity", range="A41:O46")
# VaxSeverTotal<-data.frame(VaxSeverTotal)
# VaxSeverTotal<-VaxSeverTotal[-c(2:3)]
# colnames(VaxSeverTotal)<-c("Total",	"Type", "Severe",	"Moderate"	, "Not Severe (mild)",	"Prop Severe	Lower CI",	"Upper CI",	"Percent_Severe",	"Percent_Lower_CI",	"Percent_Upper_CI",	"Study Size",	"Weighting")
# VaxSeverTotal <- VaxSeverTotal %>% filter(Total=="Weighted") %>% dplyr::select(Percent_Severe, Percent_Lower_CI, Percent_Upper_CI, Type, Total)
# 
# VaxLong<-melt(VaxSever, id.vars = c("Immunization.Status", "Location", "Study.Size"))
# VaxTotalLong<-melt(VaxSeverTotal, id.vars = c("Type", "Total"))
# 
# dodge <- position_dodge(width=0.45)
# 
# ggplot() + geom_point(data=VaxSever, aes(x=Immunization.Status, y=Percent.Severe, group=Location, size = Study.Size, color=Location), shape=15, position=dodge) + geom_errorbar(data=VaxSever, mapping=aes(x=Immunization.Status, ymin=Percent.Lower.CI, ymax=Percent.Upper.CI,color=Location), width=0.1, size=0.6, position=dodge) + geom_violin(data=VaxTotalLong, aes(x=Type, y=value, group=interaction(Type)),position=position_nudge(x=0.33), width=0.15, fill="gray30",trim = FALSE) + xlab("Immunization Status") + ylab("Severe Cases (%)")+ scale_x_discrete(labels=c("Full" = "Full\nDTP3", "Imperfect" = "Imperfect\nDTP 1 or 2","None" = "None")) +labs(size="Study\nSize") +ylim(0,100)
# 
# 
# 
