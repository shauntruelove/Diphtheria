
# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('devtools')) install.packages('devtools'); library(devtools)
if(!require('ggbeeswarm')) install.packages('ggbeeswarm'); library(ggbeeswarm)
if(!require('ggpubr')) install.packages('ggpubr'); library(rjags)
if(!require('cowplot')) install.packages('cowplot'); library(cowplot)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('boot')) install.packages('boot'); library(boot)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)
if(!require('lme4')) install.packages('lme4'); library(lme4)
if(!require('effects')) install.packages('effects'); library(effects)
if(!require('arm')) install.packages('arm'); library(arm)

#devtools::install_github("kassambara/ggpubr")

# Get rid of scientific notation numbers
options(scipen=999)


# Load & Clean Data ---------------------------------------------------------------

CFR.DAT.data <- read.csv(file="./data/cfr_antitoxin_delay.csv", header = TRUE, stringsAsFactors = FALSE)

CFR.DAT.data <- CFR.DAT.data[, c("Ref_index", "Location", "Country", "Region",
                         "Decade", "Mid.Date", "TimeToDAT","TimeMin","TimeMax", "Cases", "Deaths")]
CFR.DAT.data$CFR <- CFR.DAT.data$Deaths/CFR.DAT.data$Cases
CFR.DAT.data$TreatDays <- CFR.DAT.data$TimeMax - CFR.DAT.data$TimeMin + 1

# CFR by DAT Delay individualized

CFR.DAT.long <-  CFR.DAT.data %>% mutate(row = 1:nrow(CFR.DAT.data)) %>%
    group_by(row) %>%
    do({
        left_join(data_frame(row=.$row, died = c(rep(0,.$Cases-.$Deaths), rep(1,.$Deaths))),.,by='row')
    })

CFR.DAT.long$TreatDelay <- NULL
CFR.DAT.long$case <- 1


# Set up Delay groups
CFR.DAT.long$TreatDelay <- CFR.DAT.long$TimeMin
CFR.DAT.long$TreatDelay[which(CFR.DAT.long$TimeMin>=5)] <- "5+"
CFR.DAT.long <- CFR.DAT.long[!(CFR.DAT.long$TimeMin<5 & CFR.DAT.long$TimeMax>CFR.DAT.long$TimeMin),]

# Sum cases and deaths for the split ages back to the new groups & Calc CFR
CFR.DAT.grouped <- CFR.DAT.long %>% group_by(Ref_index, Location, Region, Decade, TreatDelay) %>%
    summarise(Cases=round(sum(case, na.rm=TRUE),0), Deaths=round(sum(died, na.rm=TRUE),0)) %>% as.data.frame
CFR.DAT.grouped$CFR <- round(CFR.DAT.grouped$Deaths / CFR.DAT.grouped$Cases,2)
CFR.DAT.grouped$CFR[CFR.DAT.grouped$Deaths==0] <- 0
CFR.DAT.grouped$CFR <- CFR.DAT.grouped$CFR * 100



# Basic Pooled Analysis - Logistic Regression - No Mixed Effect -----------

mod <- glm(died~TreatDelay, data=CFR.DAT.long, family = binomial(link = "logit"))

probs1 <- predict(mod, data.frame(TreatDelay=unique(CFR.DAT.long$TreatDelay)), type = "response")
names(probs1) <- unique(CFR.DAT.long$TreatDelay)

prob1.se <- predict(mod, data.frame(TreatDelay=unique(CFR.DAT.long$TreatDelay)), se.fit=TRUE)
probs1 <- cbind(exp(prob1.se$fit) / (1+exp(prob1.se$fit)),
                exp(prob1.se$fit - 1.96*prob1.se$se.fit) / (1+exp(prob1.se$fit - 1.96*prob1.se$se.fit)),
                exp(prob1.se$fit + 1.96*prob1.se$se.fit) / (1+exp(prob1.se$fit + 1.96*prob1.se$se.fit)))
row.names(probs1) <- unique(CFR.DAT.long$TreatDelay)



# Model with mixed effect for Study (we used this one in the paper) ---------------------------------------
CFR.DAT.long$Ref_index <- factor(CFR.DAT.long$Ref_index)
CFR.DAT.long$TreatDelay <- factor(CFR.DAT.long$TreatDelay, levels = c("1","2","3","4","5+"))

mod2 <- glmer(died ~ TreatDelay + (1|Ref_index), data=CFR.DAT.long, family = binomial(link = "logit"))
summary(mod2)
exp(fixef(mod2))
exp(cbind("Odds ratio" = fixef(mod2), confint(mod2, method="Wald", level = 0.95)[-1,]))
se.fixef (mod2)

se <- sqrt(diag(vcov(mod2)))
# table of estimates with 95% CI
(betas.df <- cbind(Est = fixef(mod2), LL = fixef(mod2) - 1.96 * se, UL = fixef(mod2) + 1.96 * se))

rbind(sapply(betas.df[1,], inv.logit),
      sapply(c(betas.df[1,1]+betas.df[2,1], betas.df[1,2]+betas.df[2,2],  betas.df[1,3]+betas.df[2,3]), inv.logit),
      sapply(c(betas.df[1,1]+betas.df[3,1], betas.df[1,2]+betas.df[3,2],  betas.df[1,3]+betas.df[3,3]), inv.logit),
      sapply(c(betas.df[1,1]+betas.df[4,1], betas.df[1,2]+betas.df[4,2],  betas.df[1,3]+betas.df[4,3]), inv.logit),
      sapply(c(betas.df[1,1]+betas.df[5,1], betas.df[1,2]+betas.df[5,2],  betas.df[1,3]+betas.df[5,3]), inv.logit))

# Convert to Probabilities (verifying the above with effects package)
(probs2 <- data.frame(effect(c("TreatDelay"), mod2)))



# Use log-binomial model to explicitly calculate the RR ---------------------------------------
mod2b <- glmer(died ~ TreatDelay + (1|Ref_index), data=CFR.DAT.long, family = poisson(link = "log"))
summary(mod2b)
exp(fixef(mod2b))
exp(cbind("Risk ratio" = fixef(mod2b), confint(mod2b, method="Wald", level = 0.95)[-1,]))[-1,]

# Convert to Probabilities
(ef.1 <- data.frame(effect(c("TreatDelay"), mod2)))
probs1



# Get Numbers of studies and subjects for Table 2 ---------------------------------------------------

length(unique(CFR.DAT.long$Ref_index))
nrow(CFR.DAT.long)





# BOOTSTRAPPED RESULTS ----------------------------------------------------
#  -- From https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/

n.boots <- 500
run.boots <- FALSE  

if (run.boots){
  
  sampler <- function(dat, clustervar, replace = TRUE, reps = 1) {
      cid <- unique(dat[, clustervar[1]])
      ncid <- length(cid)
      recid <- sample(cid, size = ncid * reps, replace = TRUE)
      if (replace) {
          rid <- lapply(seq_along(recid), function(i) {
              cbind(NewID = i, RowID = sample(which(dat[, clustervar] == recid[i]), size = length(which(dat[, clustervar] == recid[i])), replace = TRUE))
          })
      } else {
          rid <- lapply(seq_along(recid), function(i) { cbind(NewID = i, RowID = which(dat[, clustervar] == recid[i])) })
      }
      dat <- as.data.frame(do.call(rbind, rid))
      dat$Replicate <- factor(cut(dat$NewID, breaks = c(1, ncid * 1:reps), include.lowest = TRUE, labels = FALSE))
      dat$NewID <- factor(dat$NewID)
      return(dat)
  }
  
  # Using mod2 as final model
  #mod2 <- glmer(died ~ TreatDelay + (1|Ref_index), data=CFR.DAT.long, family = binomial(link = "logit"))
  dat.indiv <- CFR.DAT.long
  
  set.seed(20)
  dat.indiv$ID <- 1:nrow(dat.indiv)
  dat.indiv <- data.frame(dat.indiv[, c('Ref_index','died','TreatDelay')])
  tmp <- sampler(dat.indiv, "Ref_index", reps=n.boots)
  bigdata <- cbind(tmp, dat.indiv[tmp$RowID,])
  tmp <- NULL
  
  print(object.size(bigdata), units='MB')
  save(dat.indiv, bigdata, mod2, file = paste0('data/cfr_antitoxin_delay_bigdata_', n.boots, 'boots.RData'))
  
  f <- fixef(mod2)
  r <- getME(mod2, "theta")
  
  # Run the bootstrap procedure in parallel to save time
  library(parallel)
  cl <- makeCluster(round((parallel::detectCores()) / 2)) # Make a cluster using half of available cores
  clusterExport(cl, c("bigdata", "f", "r"))
  clusterEvalQ(cl, require(lme4))
  
  # Bootstrap function, with ME logistic regression model
  myboot <- function(i) {
      object <- try(glmer(died ~ TreatDelay + (1 | NewID), data = bigdata, subset = Replicate==i, family=binomial(link='logit'),
                          nAGQ = 1, start = list(fixef = f, theta = r)), silent = TRUE)
      if (class(object) == "try-error") {
          return(object)
      }
      c(fixef(object), getME(object, "theta"))
  }
  
  # Run the bootstrap
  start <- proc.time()
  res <- parLapplyLB(cl, X = levels(bigdata$Replicate), fun = myboot)
  end <- proc.time()
  end-start
  
  # shut down the cluster
  stopCluster(cl)
  
  # calculate proportion of models that successfully converged
  success <- sapply(res, is.numeric)
  mean(success)
  
  # combine successful results
  bigres <- do.call(cbind, res[success])
  
  save(bigres, file=paste0('results/boots_antitoxin_delay_', n.boots, 'boots.RData'))
  rm(res)
}

# Load Data if Saved from Bootstrap
#load(file='results/boots_antitoxin_delay_500boots.RData') # Loads bigres
#load(file='results/boots_antitoxin_delay.RData') # Loads bigres 100 bootstraps


# calculate 2.5th and 97.5th percentiles for 95% CI
(ci <- t(apply(bigres, 1, quantile, probs = c(0.025, 0.975))))

# All results
finaltable <- as.data.frame(cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres), ci))
# round and print
(finaltable <- round(finaltable, 3))

inv.logit(finaltable)


# PROBABILITIES -----------------------------------------------------------

boot.prob.1 <- inv.logit(bigres[1,])
boot.prob.2 <- inv.logit(bigres[1,]+bigres[2,])
boot.prob.3 <- inv.logit(bigres[1,]+bigres[3,])
boot.prob.4 <- inv.logit(bigres[1,]+bigres[4,])
boot.prob.5 <- inv.logit(bigres[1,]+bigres[5,])


probs_boot <- rbind(c(mean(boot.prob.1), quantile(boot.prob.1, probs=c(.5,.025,.975,.25,.75)), sd(boot.prob.1) ),
                    c(mean(boot.prob.2), quantile(boot.prob.2, probs=c(.5,.025,.975,.25,.75)), sd(boot.prob.2) ),
                    c(mean(boot.prob.3), quantile(boot.prob.3, probs=c(.5,.025,.975,.25,.75)), sd(boot.prob.3) ),
                    c(mean(boot.prob.4), quantile(boot.prob.4, probs=c(.5,.025,.975,.25,.75)), sd(boot.prob.4) ),
                    c(mean(boot.prob.5), quantile(boot.prob.5, probs=c(.5,.025,.975,.25,.75)), sd(boot.prob.5) ))
colnames(probs_boot) <- c('mean.prob','median.prob','LL','UL','25pct','75pct','SD')
row.names(probs_boot) <- c(1,2,3,4,"5+")
probs_boot


write.csv(probs_boot, file='results/cfr_antitoxin_delay_boots_results_table.csv', row.names = FALSE)

probability.boots.raw <- rbind(data.frame(prob=boot.prob.1, delay='1'),
                               data.frame(prob=boot.prob.2, delay='2'),
                               data.frame(prob=boot.prob.3, delay='3'),
                               data.frame(prob=boot.prob.4, delay='4'),
                               data.frame(prob=boot.prob.5, delay='5+'))

write.csv(probability.boots.raw, file='results/cfr_antitoxin_delay_boots_results_raw.csv', row.names = FALSE)





