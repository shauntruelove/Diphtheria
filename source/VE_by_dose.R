
# CFR Estimated using Stan ------------------------------------------------


# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('caret')) install.packages('caret'); library(caret)
if(!require('rstan')) install.packages('rstan'); library(rstan)
library(ggplot2)

source(file.path("source/VE_data_source.R"))

# Get rid of scientific notation numbers
options(scipen=999)

#......................................................................................
# VE by vacc --------------------------------------------------------------
#......................................................................................

VE_long <- load_VE_long() # Run load script to load and clean VE data
VE_long$Dose013 <- as.character(VE_long$Dose013)
VE_long_full <- VE_long %>% filter(Dose013=='None' | Dose013=='Full')
VE_long_partial <- VE_long %>% filter(Dose013=='None' | Dose013=='Partial')
VE_long_full$Dose013 <- relevel(as.factor(VE_long_full$Dose013), ref='None')
VE_long_partial$Dose013 <- relevel(as.factor(VE_long_partial$Dose013), ref='None')

dmy <- dummyVars(case ~ Dose013, data=VE_long_full)
fixed.dummies <- data.frame(predict(dmy, newdata = VE_long_full))
dmy <- dummyVars(case ~ Dose013, data=VE_long_partial)
fixed.dummies.partial <- data.frame(predict(dmy, newdata = VE_long_partial))


VE_long_partial$Study <- as.integer(as.factor(VE_long_partial$Study))
VE_long_full$Study <- as.integer(as.factor(VE_long_full$Study))

VE.full.stan <- list(N=as.integer(nrow(VE_long_full)),
                     J=as.integer(length(unique(VE_long_full$Study))),
                     id=as.integer(VE_long_full$Study),
                     K=2,
                     X=as.matrix(fixed.dummies),
                     Y=as.integer(VE_long_full$case))

VE.partial.stan <- list(N=as.integer(nrow(VE_long_partial)),
                     J=as.integer(length(unique(VE_long_partial$Study))),
                     id=as.integer(VE_long_partial$Study),
                     K=2,
                     X=as.matrix(fixed.dummies.partial),
                     Y=as.integer(VE_long_partial$case))


# MODEL ---------------------------------------------------------------

fileName <- "./source/VE_mixed_logistic.stan"
ret <- stanc(fileName) # Check Stan file
ret_sm <- stan_model(stanc_ret = ret) # Compile Stan code
save(fileName, ret, ret_sm, file="./source/VE_mixed_logistic_compiled.RData")
load(file="./source/VE_mixed_logistic_compiled.RData")


# RUN TEST ---------------------------------------------------------------
rstan_options (auto_write=TRUE)
options (mc.cores=1) # Run on multiple cores

test_fit <- sampling(ret_sm, warmup=50, iter=100, seed=123, data=VE.full.stan, chains=1, control=list(adapt_delta=.85))
test_fit
summary(rstan::extract(test_fit)$VE)
quantile(rstan::extract(test_fit)$VE, probs = c(0.025,0.975))


# Run Model - Full Vaccination ---------------------------------------------------------------
rstan_options (auto_write=TRUE)
options (mc.cores=3) # Run on multiple cores

VE.full.stan.fit <- sampling(ret_sm, warmup=1000, iter=2000, seed=123, data=VE.full.stan, 
                        chains=3, control=list(adapt_delta=.95))
VE.full.stan.fit

save(VE.full.stan.fit, file='./results/VE.full.stan.RData')
load(file='./results/VE.full.stan.RData') # VE.full.stan.fit

summary(rstan::extract(VE.full.stan.fit)$VE)
quantile(rstan::extract(VE.full.stan.fit)$VE, probs = c(0.025,0.975))

chains <- rstan::extract(VE.full.stan.fit)

VE_all <- 1- exp(chains$beta)
summary(VE_all)
quantile(VE_all, probs=c(0.025,.975))

VE_study <- 1- exp(chains$beta_j)
summary(VE_study)
    



    
# Run Model - Partial Vaccination ---------------------------------------------------------------
rstan_options (auto_write=TRUE)
options (mc.cores=3) # Run on multiple cores

VE.partial.stan.fit <- sampling(ret_sm, warmup=500, iter=1000, seed=123, data=VE.partial.stan, 
                             chains=3, control=list(adapt_delta=.95))
VE.partial.stan.fit

save(VE.partial.stan.fit, file='./results/VE.partial.stan.RData')
load(file='./results/VE.partial.stan.RData') # VE.partial.stan.fit

summary(rstan::extract(VE.partial.stan.fit)$VE)
quantile(rstan::extract(VE.partial.stan.fit)$VE, probs = c(0.025,0.975))

chains=extract(rstan::VE.partial.stan.fit)

VE_all <- 1- exp(chains$beta)
summary(VE_all)
quantile(VE_all, probs=c(0.025,.975))

VE_study <- 1- exp(chains$beta_j)
summary(VE_study)
tmp <- apply(chains$beta_j, 2, quantile, probs=c(0.025,0.975))
1-exp(tmp)







# Check with metafor package ----------------------------------------------


VE_full <- VE_long_full[!duplicated(VE_long_full$row),]
VE_full <- VE_full %>% group_by(Population, Dose013) %>% summarise(Cases=sum(Cases_total), Controls=sum(Controls))
VE_full_wide <- VE_full %>% select(-Controls) %>% spread(key=Dose013, value = Cases) %>% rename(Case_None=None, Case_Full=Full)
VE_full_wide2 <- VE_full %>% select(-Cases) %>% spread(key=Dose013, value = Controls) %>% rename(Control_None=None, Control_Full=Full)
VE_full_wide <- full_join(VE_full_wide, VE_full_wide2)

library(metafor)
dat <- escalc(measure = "OR", ai = Case_Full, bi = Control_Full, 
              ci = Case_None, di = Control_None, 
              data = VE_full_wide, append = TRUE)
res <- rma(yi, vi, data = dat)
res


VE_partial <- VE_long_partial[!duplicated(VE_long_partial$row),]
VE_partial <- VE_partial %>% group_by(Population, Dose013) %>% summarise(Cases=sum(Cases_total), Controls=sum(Controls))
VE_partial_wide <- VE_partial %>% select(-Controls) %>% spread(key=Dose013, value = Cases) %>% rename(Case_None=None, Case_Partial=Partial)
VE_partial_wide2 <- VE_partial %>% select(-Cases) %>% spread(key=Dose013, value = Controls) %>% rename(Control_None=None, Control_Partial=Partial)
VE_partial_wide <- full_join(VE_partial_wide, VE_partial_wide2)

library(metafor)
dat <- escalc(measure = "OR", ai = Case_Partial, bi = Control_Partial, 
              ci = Case_None, di = Control_None, 
              data = VE_partial_wide, append = TRUE)
res <- rma(yi, vi, data = dat)
res


