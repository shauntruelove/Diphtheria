
# VE by Age ------------------------------------------------

# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)

# Get rid of scientific notation numbers
options(scipen=999)



#......................................................................................
# VE by vacc --------------------------------------------------------------
#......................................................................................

VE_long <- read.csv(file="./data/VE_age.csv", header = TRUE, stringsAsFactors = FALSE)
VE_long$Dose013 <- factor(VE_long$Dose013, levels=c('None','Partial','Full'))

VE_long0_4 <- VE_long %>% filter(AgeGroup3=='0-4y')
VE_long5_19 <- VE_long %>% filter(AgeGroup3=='5-19y')
VE_long20 <- VE_long %>% filter(AgeGroup3=='20+')


# Calc using glm ---------------------------------------------------
# - because age groups are separately available from different sources, separate regressions will be done.

# 0-4y
lm1 <- glm(case ~ Dose013, data=VE_long0_4, family=binomial(link='logit'))
coefs1 <- coef(lm1)
1-exp(coefs1[2:3])
1-exp(confint(lm1))

# 5-19y
lm2 <- glm(case ~ Dose013, data=VE_long5_19, family=binomial(link='logit'))
coefs2 <- coef(lm2)
1-exp(coefs2[2:3])
1-exp(confint(lm2))

# >=20y
lm3 <- glm(case ~ Dose013, data=VE_long20, family=binomial(link='logit'))
coefs3 <- coef(lm3)
1-exp(coefs3[2:3])
1-exp(confint(lm3))

