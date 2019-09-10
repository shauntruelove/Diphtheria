

# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('devtools')) install.packages('devtools'); library(devtools)
if(!require('ggpubr')) install.packages('ggpubr'); library(rjags)
if(!require('cowplot')) install.packages('cowplot'); library(cowplot)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('boot')) install.packages('boot'); library(boot)
if(!require('caret')) install.packages('caret'); library(caret)
library(brms)
library(ggplot2)
library(gridExtra)
library(lme4)

# Get rid of scientific notation numbers
options(scipen=999)



#......................................................................................
# CFR by vaccination ------------------------------------------------------
#......................................................................................

CFR.Vax.data.orig <- read.csv(file='./data/cfr_vaccination.csv', header = TRUE, stringsAsFactors = FALSE)
CFR.Vax.data.orig <- CFR.Vax.data.orig[!is.na(CFR.Vax.data.orig$Case.none) & !is.na(CFR.Vax.data.orig$Study),] # Remove Studies with missing data

CFR.Vax.data <- CFR.Vax.data.orig
CFR.Vax.data$CFR.full <- CFR.Vax.data$CFR.full*100
CFR.Vax.data$CFR.partial <- CFR.Vax.data$CFR.partial*100
CFR.Vax.data$CFR.non <- CFR.Vax.data$CFR.non*100

# Melt to long form
CFR.Vax.data.tmp <- CFR.Vax.data[,c("Ref_index", "Location", "Country", "Region", "Decade", "Year.Mid", 
                                    "CFR.full", "CFR.partial","CFR.non")]
CFR.Vax <- melt(CFR.Vax.data.tmp, id.vars = c("Ref_index", "Location", "Country", "Region", "Decade", "Year.Mid"))
CFR.Vax <- CFR.Vax[!is.na(CFR.Vax$value),]

# Mean and 95% CI
CFR.Vax %>% filter(variable=="CFR.full") %>% summarize(mean=mean(value), ql = quantile(value,.025),qu = quantile(value,.975))
CFR.Vax %>% filter(variable=="CFR.partial") %>% summarize(mean=mean(value), ql = quantile(value,.025),qu = quantile(value,.975))
CFR.Vax %>% filter(variable=="CFR.non") %>% summarize(mean=mean(value), ql = quantile(value,.025),qu = quantile(value,.975))

#Median & IQR
CFR.Vax %>% filter(variable=="CFR.full") %>% summarize(median=median(value), ql = quantile(value,.25),qu = quantile(value,.75))
CFR.Vax %>% filter(variable=="CFR.partial") %>% summarize(median=median(value), ql = quantile(value,.25),qu = quantile(value,.75))
CFR.Vax %>% filter(variable=="CFR.non") %>% summarize(median=median(value), ql = quantile(value,.25),qu = quantile(value,.75))

p3 <- ggplot(CFR.Vax, aes(x=variable, y=value)) +geom_boxplot() +xlab("Vaccination Status")+ylab("Case Fatality Ratio") + scale_x_discrete(labels=c("None" = "None", "Partial" = "Imperfect\nDTP1 or 2", "Full" = "Full\nDTP3"))
p3



# CFR by vaccination OR -----------------------------------------------------------

CFR.Vax <- CFR.Vax.data.orig
CFR.Vax.1 <- CFR.Vax[,c("Study", "Region", "Decade", "Case.full", "Case.partial", "Case.none")]
colnames(CFR.Vax.1) <- c("Study", "Region", "Decade", "Full", "Partial", "None")
CFR.Vax.2 <- CFR.Vax[,c("Study", "Region", "Decade", "Death.full", "Death.partial", "Death.none")]
colnames(CFR.Vax.2) <- c("Study", "Region", "Decade", "Full", "Partial", "None")

CFR.Vax.long1 <- melt(CFR.Vax.1, id.vars = c("Study", "Decade", "Region"))
colnames(CFR.Vax.long1) <- c("Study", "Region", "Decade", "variable", "Cases")
CFR.Vax.long2 <- melt(CFR.Vax.2, id.vars = c("Study", "Decade", "Region"))
colnames(CFR.Vax.long2) <- c("Study", "Region", "Decade", "variable", "Deaths")

CFR.Vax.long <- cbind(CFR.Vax.long1, CFR.Vax.long2$Deaths)
colnames(CFR.Vax.long) <- c("Study", "Region", "Decade", "Vacc","Cases" ,"Deaths")
rm(CFR.Vax.1, CFR.Vax.2, CFR.Vax.long1, CFR.Vax.long2)

CFR.Vax.long <- CFR.Vax.long[CFR.Vax.long$Cases!=0,]
CFR.Vax.long <- CFR.Vax.long %>% mutate(row = 1:nrow(CFR.Vax.long))
CFR.Vax.long <- CFR.Vax.long %>%
  group_by(row) %>%
  do({
    left_join(data_frame(row=.$row,died = c(rep(0,.$Cases-.$Deaths),rep(1,.$Deaths))),.,by='row')
  })



mod <- glm(died ~ Vacc, data=CFR.Vax.long, family = binomial(link = "logit"))
exp(coef(mod))[2]
exp(cbind("Odds ratio" = coef(mod), confint.default(mod, level = 0.95)))

# Calculate the probabilities
probs1 <- predict(mod, data.frame(Vacc=unique(CFR.Vax.long$Vacc)), type = "response")
names(probs1) <- unique(CFR.Vax.long$Vacc)
probs1

prob1.se <- predict(mod, data.frame(Vacc=unique(CFR.Vax.long$Vacc)), se.fit=TRUE)

probs1 <- cbind(exp(prob1.se$fit) / (1+exp(prob1.se$fit)),
                exp(prob1.se$fit - 1.96*prob1.se$se.fit) / (1+exp(prob1.se$fit - 1.96*prob1.se$se.fit)),
                exp(prob1.se$fit + 1.96*prob1.se$se.fit) / (1+exp(prob1.se$fit + 1.96*prob1.se$se.fit)))
row.names(probs1) <- unique(CFR.Vax.long$Vacc)
probs1





# Model with mixed effect for Study  ---------------------------------------
CFR.Vax.long$Study <- factor(CFR.Vax.long$Study)

mod3 <- glmer(died ~ Vacc + (1|Study), data=CFR.Vax.long, family = binomial(link = "logit"))
summary(mod3)
exp(fixef(mod3))
exp(cbind("Odds ratio" = fixef(mod3), confint(mod3, method="Wald", level = 0.95)[-1,]))

# Convert to Probabilities
library(effects)
(ef.1 <- data.frame(effect(c("Vacc"), mod3)))





#........................................................................................
# Vaccination & CFR - Stan Estimation -------------------------------------------------------
#........................................................................................

CFR.Vax.long_test <- CFR.Vax.long[as.integer(CFR.Vax.long$Study)>=1 & as.integer(CFR.Vax.long$Study)<=3, ]  # subset to 3 studies
dmy <- dummyVars(died ~ Vacc, data=CFR.Vax.long_test)
fixed.dummies_test <- data.frame(predict(dmy, newdata=CFR.Vax.long_test))

dmy <- dummyVars(died ~ Vacc, data=CFR.Vax.long)
fixed.dummies <- data.frame(predict(dmy, newdata = CFR.Vax.long))

CFR.Vax.stan_test <- list(N=as.integer(nrow(CFR.Vax.long_test)),
                          J=as.integer(length(unique(CFR.Vax.long_test$Study))),
                          id=as.integer(CFR.Vax.long_test$Study),
                          K=3,
                          X=as.matrix(fixed.dummies_test),
                          Y=as.integer(CFR.Vax.long_test$died))

CFR.Vax.stan <- list(N=as.integer(nrow(CFR.Vax.long)),
                     J=as.integer(length(unique(CFR.Vax.long$Study))),
                     id=as.integer(CFR.Vax.long$Study),
                     K=3,
                     X=as.matrix(fixed.dummies),
                     Y=as.integer(CFR.Vax.long$died))

table(CFR.Vax.long$Study)


# MODEL ---------------------------------------------------------------

fileName <- "source/CFR_mixed_logistic.stan"
ret <- stanc(fileName) # Check Stan file
ret_sm <- stan_model(stanc_ret = ret) # Compile Stan code
save(fileName, ret, ret_sm, file="source/CFR_mixed_logistic_compiled.Rdata")
load(file="source/CFR_mixed_logistic_compiled.Rdata")



# RUN TEST ---------------------------------------------------------------
rstan_options (auto_write=TRUE)
options (mc.cores=1) # Run on multiple cores

start_time <- Sys.time()
test_fit <- sampling(ret_sm, warmup=50, iter=100, seed=123, data=CFR.Vax.stan_test, 
                     chains=1, control=list(adapt_delta=.85))
end_time <- Sys.time()
end_time - start_time
test_fit

summary(extract(test_fit)$CFR0)
summary(extract(test_fit)$CFR1)
summary(extract(test_fit)$CFR2)



# Run Model ---------------------------------------------------------------
rstan_options (auto_write=TRUE)
options (mc.cores=3) # Run on multiple cores

CFR.Vax.stan.fit <- sampling(ret_sm, warmup=500, iter=1000, seed=123, data=CFR.Vax.stan, 
                             chains=3, control=list(adapt_delta=.95))

save(CFR.Vax.stan.fit, file='results/CFR.Vax.stan.Study.Rdata')
load(file='results/CFR.Vax.stan.Study.Rdata') # CFR.Vax.stan.fit

CFR.Vax.stan.fit

stan.fits <- rstan::extract(CFR.Vax.stan.fit)

CFR0 <- stan.fits$CFR0 # No vaccination
CFR1 <- stan.fits$CFR1 # Full Vaccination
CFR2 <- stan.fits$CFR2 # Partial Vaccination

cfr_vax <- data.frame(NoVacc=CFR0, Partial=CFR2, Full=CFR1)
write.csv(cfr_vax, file='results/cfr_vax.csv', row.names = FALSE)

# Unvaccinated
summary(CFR0); quantile(CFR0, probs = c(0.025, 0.975))
# Partially Vaccinated
summary(CFR2); quantile(CFR2, probs = c(0.025, 0.975))
# Fully Vaccinated
summary(CFR1); quantile(CFR1, probs = c(0.025, 0.975))


# Relative Risks
RR0to1 <-  CFR1 / CFR0
summary(RR0to1); quantile(RR0to1, probs = c(0.025, 0.975))

RR0to2 <-  CFR2 / CFR0
summary(RR0to2); quantile(RR0to2, probs = c(0.025, 0.975))








