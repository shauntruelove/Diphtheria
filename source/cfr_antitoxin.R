
# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('caret')) install.packages('caret'); library(caret)
if(!require('rstan')) install.packages('rstan'); library(rstan)
if(!require('brms')) install.packages('brms'); library(brms)
if(!require('ggbeeswarm')) install.packages('ggbeeswarm'); library(ggbeeswarm)
if(!require('lme4')) install.packages('lme4'); library(lme4)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)


# Get rid of scientific notation numbers
options(scipen=999)


#........................................................................................
# ANTITOXIN RECEIPT -------------------------------------------------------
#........................................................................................

# LOAD AND CLEAN DATA -----------------------------------------------------

CFR.DAT <- read.csv(file="data/cfr_antitoxin.csv", header = TRUE, stringsAsFactors = FALSE)
CFR.DAT <- CFR.DAT[!is.na(CFR.DAT$Cases),]
CFR.DAT$CFR <- (CFR.DAT$Deaths/CFR.DAT$Cases)*100
CFR.DAT <- CFR.DAT[!is.na(CFR.DAT$Cases),]
CFR.DAT$LocDecade <- paste0(CFR.DAT$Location, ' - ', CFR.DAT$Decade)

# Covert to long
CFR.DAT.long <- CFR.DAT %>% mutate(row = 1:nrow(CFR.DAT))
CFR.DAT.long <- CFR.DAT.long %>% group_by(row) %>%
  do({ left_join(data_frame(row=.$row,died = c(rep(0,.$Cases-.$Deaths),rep(1,.$Deaths))),.,by='row') })

CFR.DAT.long$Antitoxin <- factor(CFR.DAT.long$Antitoxin)
CFR.DAT.long$Decade <- factor(CFR.DAT.long$Decade)
CFR.DAT.long$Location <- factor(CFR.DAT.long$Location)
CFR.DAT.long$LocDecade <- factor(CFR.DAT.long$LocDecade)
CFR.DAT.long$LocStudy <- paste0(CFR.DAT.long$Location,"-",CFR.DAT.long$Ref_index)
CFR.DAT.long$Ref_index[is.na(CFR.DAT.long$Ref_index)] <- paste0(CFR.DAT.long$Location[is.na(CFR.DAT.long$Ref_index)],"-",CFR.DAT.long$Decade[is.na(CFR.DAT.long$Ref_index)])
CFR.DAT.long$Ref_index <- as.factor(CFR.DAT.long$Ref_index)
CFR.DAT.long$LocStudy <- as.factor(CFR.DAT.long$LocStudy)




# BASIC ANALYSIS ----------------------------------------------------------

mod <- glm(died~Antitoxin, data=CFR.DAT.long, family = poisson(link = "log"))
exp(coef(mod))[2]
exp(cbind("Risk ratio" = coef(mod), confint.default(mod, level = 0.95)))

p <- ggplot(CFR.DAT, aes(x=Antitoxin, y=CFR)) +geom_boxplot()+ 
  geom_beeswarm(data=CFR.DAT, aes(color=Decade, shape=Region), size=2) + 
  xlab("Antitoxin")+ylab("Case Fatality Ratio") +xlab("Antitoxin")+ylab("Case Fatality Ratio") +
  #scale_shape_manual(values=c(16, 24))+scale_colour_gradient( limits=c(1730, 2010)) + 
  guides(shape = guide_legend(order=1),colour = guide_colourbar(order = 2))  
p


# Calculate the probabilities
probs1 <- predict(mod, data.frame(Antitoxin=unique(CFR.DAT.long$Antitoxin)), type = "response")
names(probs1) <- unique(CFR.DAT.long$Antitoxin)
probs1

prob1.se <- predict(mod, data.frame(Antitoxin=unique(CFR.DAT.long$Antitoxin)), se.fit=TRUE)
probs1 <- cbind(exp(prob1.se$fit) / (1+exp(prob1.se$fit)),
                exp(prob1.se$fit - 1.96*prob1.se$se.fit) / (1+exp(prob1.se$fit - 1.96*prob1.se$se.fit)),
                exp(prob1.se$fit + 1.96*prob1.se$se.fit) / (1+exp(prob1.se$fit + 1.96*prob1.se$se.fit)))
row.names(probs1) <- unique(CFR.DAT.long$Antitoxin)
probs1



#........................................................................................
# ANTITOXIN RECEIPT - Stan Estimation -------------------------------------------------------
#........................................................................................

CFR.DAT.long_test <- CFR.DAT.long[as.integer(CFR.DAT.long$LocDecade)>=1 & as.integer(CFR.DAT.long$LocDecade)<=3, ]  # subset to 3 studies
dmy <- dummyVars(died ~ Antitoxin, data=CFR.DAT.long_test)
fixed.dummies_test <- data.frame(predict(dmy, newdata=CFR.DAT.long_test))

dmy <- dummyVars(died ~ Antitoxin, data=CFR.DAT.long)
fixed.dummies <- data.frame(predict(dmy, newdata = CFR.DAT.long))

CFR.DAT.stan_test <- list(N=as.integer(nrow(CFR.DAT.long_test)),
                          J=as.integer(length(unique(CFR.DAT.long_test$LocStudy))),
                          id=as.integer(CFR.DAT.long_test$LocStudy),
                          K=2,
                          X=as.matrix(fixed.dummies_test),
                          Y=as.integer(CFR.DAT.long_test$died))

CFR.DAT.stan <- list(N=as.integer(nrow(CFR.DAT.long)),
                     J=as.integer(length(unique(CFR.DAT.long$LocStudy))),
                     id=as.integer(CFR.DAT.long$LocStudy),
                     K=2,
                     X=as.matrix(fixed.dummies),
                     Y=as.integer(CFR.DAT.long$died))

table(CFR.DAT.long$LocStudy)
table(CFR.DAT.long$Decade)


# MODEL ---------------------------------------------------------------

fileName <- "source/CFR_mixed_2lvl_logistic.stan"
ret <- stanc(fileName) # Check Stan file
ret_sm <- stan_model(stanc_ret = ret) # Compile Stan code
save(fileName, ret, ret_sm, file="source/CFR_mixed_2lvl_logistic_compiled.RData")
load(file="source/CFR_mixed_2lvl_logistic_compiled.RData")



# RUN TEST ---------------------------------------------------------------
rstan_options (auto_write=TRUE)
options (mc.cores=1) # Run on multiple cores

start_time <- Sys.time()
test_fit <- sampling(ret_sm, warmup=50, iter=100, seed=123, data=CFR.DAT.stan_test, 
                     chains=1, control=list(adapt_delta=.85))
end_time <- Sys.time()
end_time - start_time
test_fit

summary(extract(test_fit)$CFR0)
summary(extract(test_fit)$CFR1)





# Run Model ---------------------------------------------------------------
rstan_options(auto_write=TRUE)
options(mc.cores=3) # Run on multiple cores

CFR.DAT.stan.fit <- sampling(ret_sm, warmup=500, iter=1000, seed=123, data=CFR.DAT.stan, 
                             chains=3, control=list(adapt_delta=.95))

#save(CFR.DAT.stan.fit, file='results/CFR.DAT.stan.LocStudy.RData')
load(file='results/CFR.DAT.stan.LocStudy.RData') # CFR.DAT.stan.fit
#load(file='results/CFR.DAT.stan.LocDecade.RData') # CFR.DAT.stan.fit

CFR.DAT.stan.fit


CFR0 <- extract(CFR.DAT.stan.fit)$CFR0 # No antitoxin
CFR1 <- extract(CFR.DAT.stan.fit)$CFR1 # Antitoxin

cfr.antitoxin <- data.frame(antitoxin=CFR1, no_antitoxin=CFR0)
write.csv(cfr.antitoxin, file='results/cfr_antitoxin.csv', row.names = FALSE)


summary(CFR0); quantile(CFR0, probs = c(0.025, 0.975))
summary(CFR1); quantile(CFR1, probs = c(0.025, 0.975))



# Relative Risks
RR1to0 <-  CFR0 / CFR1
summary(RR1to0); quantile(RR1to0, probs = c(0.025, 0.975))

RR0to1 <-  CFR1 / CFR0
summary(RR0to1); quantile(RR0to1, probs = c(0.025, 0.975))





    


# DOUBLE-CHECK WITH BASIC CALCS -------------------------------------------

CFR.DAT.long %>% group_by(Antitoxin) %>% summarise(sum(died) / n())
tmp <- CFR.DAT.long %>% group_by(Antitoxin, Ref_index) %>% summarise(cfr=sum(died) / n())
View(tmp)
tmp <- CFR.DAT.long %>% group_by(Antitoxin, LocDecade) %>% summarise(cfr=sum(died) / n())
View(tmp)

tmp %>% group_by(Antitoxin) %>% summarise(mean(cfr))






