
# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)
if(!require('knitr')) install.packages('knitr'); library(knitr)

if(!require('googlesheets')) install.packages('googlesheets'); library(googlesheets)
if(!require('fitdistrplus')) install.packages('fitdistrplus'); library(fitdistrplus)
if(!require('coda')) install.packages('coda'); library(coda)
if(!require('rjags')) install.packages('rjags'); library(rjags)

library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(knitr)

recalc <- F

# Check if Results Exist -------------------------------------------

if (!file.exists("./trunk/results/full.fit.jags.LWW.CARRIERS.RData") | recalc) {
  
  # Set up JAGS global parameters
  n.adapt <- 10000
  iters <- 1000000
  thin <- 500
  n.chains <- 3
  
  
  # FUNCTIONS FOR JAGS FITTING ----------------------------------------------
  
  # Function for setting up JAGS
  carriage.jags.setup.fn <- function(data) {
    
    require(rjags) #installed JAGS must be 4.1 or higher
    
    #assuming observations with no upper limit have NAs in the R slot. Currently only in carr.max column
    data$carr.max[which(is.infinite(data$carr.max))] <-  NA
    
    #missing observation are just unbounded, replace with 0
    data$carr.min[which(is.na(data$carr.min))] <-  0
    
    #make a data object for JAGS
    jags.data <- list(
      CL=data$carr.min,
      CR=data$carr.max)
    
    #Let jags know the censoring situation
    #1 means interval censored
    #2 means event occured after known time
    jags.data$VSPisCensored=rep(1, length(jags.data$CL))
    jags.data$VSPisCensored[which(is.na(jags.data$CR))]=2
    
    #define variable to hold the length of the time to event for bacterial clearance
    jags.data$Y_V <- rep(NA, nrow(data))
    
    #set the intitial values for time to event (i.e., Y_V)
    VSPyInit =  jags.data$CL
    VSPyInit[which(is.na(jags.data$CL)==T)]=0
    VSPyInit[which(VSPyInit==0)]=0.000000001
    
    
    #identify which observations are interval censored and which are right censored. 
    jags.data$Vdic=which(jags.data$VSPisCensored==1)
    jags.data$Vrc=which(jags.data$VSPisCensored==2)
    
    # add to right bound to fix problem with having left and right bounds the same
    jags.data$CR <- jags.data$CR + 0.1
    
    #set the parameters we want to track
    #  - v_v     : weibull shape
    #  - scale_v : weibull scale
    parameters <- c("v_v","scale_v") #, "E")  
    
    set.seed(12345) #if this is not included, multiple initializations maybe needed.
    #initialization function for jags
    jags.inits <-  function() {
      rc <-list(
        #E = rep(0.0000000011, n.subjects), #start with a fixed E to avoid bad starting points
        v_v = runif(1,1,10),
        scale_v = runif(1,0,1), 
        Y_V=VSPyInit)
      print(rc)
      return(rc)
    }
    
    return(list(jags.data=jags.data, jags.inits=jags.inits, parameters=parameters))
  }
  
  # Function to fit JAGS model, clean, and save
  carriage.jags.fit.fn <- function(jagsfit.LWW=jagsfit.LWW, jags.setup=carr.jags.setup, group='CARRIERS',
                                   iters, thin, n.chains) {
    
    full.fit.LWW <- coda.samples(model=jagsfit.LWW, variable.names=jags.setup$parameters, 
                                 n.iter=iters, thin=thin, n.chains=n.chains)
    
    #make all of the chains a single matrix with a burnin removed
    ABC1=as.matrix(full.fit.LWW[[1]][,])
    ABC2=as.matrix(full.fit.LWW[[2]][,])
    ABC3=as.matrix(full.fit.LWW[[3]][,])
    
    ABC1=ABC1[(thin+1):(iters/thin),]
    ABC2=ABC2[(thin+1):(iters/thin),]
    ABC3=ABC3[(thin+1):(iters/thin),]
    
    #recreate MCMC object for diagnostics
    jags.fit<-list(as.mcmc(ABC1), as.mcmc(ABC2), as.mcmc(ABC3)) 
    jags.chains <- rbind(ABC1,ABC2,ABC3)
    colnames(jags.chains) <- varnames(jags.fit[[1]])
    jags.chains <- as.data.frame(jags.chains)
    jags.params <- list(iters=iters, thin=thin, n.chains=n.chains)
    jags.fit.list <- list(jags.fit=jags.fit, jags.chains=jags.chains, 
                          jags.setup=jags.setup, jags.params=jags.setup)
    save(jags.fit.list, file=paste0("./trunk/results/CarriageDuration/jags.fits.",group,".RData"))
    
    return(jags.fit.list)
  }
  
  # Download and Save Google Sheet ------------------------------------------
  diph.sheet <- gs_title("TheDiphSheet")
  
  # Duration Colonized - Weaver 1921  -------------------------------------------
  dc.data <- gs_read(ss=diph.sheet, ws="DurationClearance")
  dc.data$Treatment <- tolower(dc.data$Treatment)
  
  # Convert to long
  dc.long <- data.frame(Case=rep(dc.data$Case, times=dc.data$Num.convert),
                        Treatment=rep(dc.data$Treatment, times=dc.data$Num.convert),
                        Antitoxin=rep(dc.data$Antitoxin, times=dc.data$Num.convert),
                        carr.min=rep(dc.data$carr.min, times=dc.data$Num.convert), 
                        carr.max=rep(dc.data$carr.max, times=dc.data$Num.convert),
                        study=rep(dc.data$Study, times=dc.data$Num.convert))
  
  
  # Subset the data  -----------------------------------------------------------
  
  table(dc.data$Treatment)
  dc.case.all <- dc.long[(dc.long$Case=='Case'),]
  dc.atb <- dc.case.all[(dc.case.all$Treatment!='none'),]
  dc.pen <- dc.case.all[(dc.case.all$Treatment=='penicillin'),]
  dc.eryth <- dc.case.all[(dc.case.all$Treatment=='erythromycin'),]
  dc.none <- dc.case.all[(dc.case.all$Treatment=='none'),]
  
  table(dc.long$Case)
  dc.carr.notreat <- dc.long[(dc.long$Case=='Carrier' & dc.long$Treatment=='none'),]
  dc.case.notreat <- dc.long[(dc.long$Case=='Case' & dc.long$Treatment=='none'),]
  
  table(dc.carr.notreat$study)
  table(dc.case.notreat$study)
  table(dc.atb$study)
  
  
  #///////////////////////////////////////////////////////////////////////////////////////
  # Bayesian MCMC -----------------------------------------------------------
  #///////////////////////////////////////////////////////////////////////////////////////
  
  
  #///////////////////////////////////////////////////////////////////////////////////////
  # Carriers ----------------------------------------------------------------
  
  # Setup data and initials
  carr.jags.setup <- carriage.jags.setup.fn(dc.carr.notreat)
  
  #initialize JAGS model
  jagsfit.LWW <- jags.model(file='./trunk/source/DistributionFitCarriage.jags', 
                            data=carr.jags.setup$jags.data, 
                            inits=carr.jags.setup$jags.inits, 
                            n.chains=n.chains, quiet=F, 
                            n.adapt=n.adapt)
  
  # Fit JAGS model, clean, and save
  carr.jags.fit <- carriage.jags.fit.fn(jagsfit.LWW = jagsfit.LWW, jags.setup=carr.jags.setup, 
                                        group='CARRIERS', iters, thin, n.chains)
  
  
  #///////////////////////////////////////////////////////////////////////////////////////
  # Cases - No Antibiotics ----------------------------------------------------------------
  
  # Set up JAGS global parameters
  n.adapt <- 1000
  iters <- 100000
  thin <- 100
  n.chains <- 3
  
  case.notreat.jags.setup <- carriage.jags.setup.fn(dc.case.notreat)
  
  #initialize JAGS model
  jagsfit.LWW <- jags.model(file='./trunk/source/DistributionFitCarriage.jags', 
                            data=case.notreat.jags.setup$jags.data, 
                            inits=case.notreat.jags.setup$jags.inits, 
                            n.chains=n.chains, quiet=F, 
                            n.adapt=n.adapt)
  
  # Fit JAGS model, clean, and save
  case.notreat.jags.fit <- carriage.jags.fit.fn(jagsfit.LWW=jagsfit.LWW, jags.setup=case.notreat.jags.setup,
                                                group='CASES_NOTREAT', iters, thin, n.chains)
  
  
  #///////////////////////////////////////////////////////////////////////////////////////
  # Cases - Antibiotics ----------------------------------------------------------------
  
  case.atb.jags.setup <- carriage.jags.setup.fn(dc.atb)
  
  
  #initialize JAGS model
  jags.model.file <- ifelse(length(case.atb.jags.setup$jags.data$Vrc)==0, 
                            './trunk/source/DistributionFitCarriage_NoRC.jags', './trunk/source/DistributionFitCarriage.jags')
  jagsfit.LWW <- jags.model(file=jags.model.file, 
                            data=case.atb.jags.setup$jags.data, 
                            inits=case.atb.jags.setup$jags.inits, 
                            n.chains=n.chains, quiet=F, 
                            n.adapt=n.adapt)
  
  # Fit JAGS model, clean, and save
  case.atb.jags.fit <- carriage.jags.fit.fn(jagsfit.LWW=jagsfit.LWW, jags.setup=case.atb.jags.setup,
                                            group='CASES_ATB', iters, thin, n.chains)
  
  
}



# JAGS Results ------------------------------------------------------------

load("./trunk/results/CarriageDuration/jags.fits.CARRIERS.RData")
carr.jags.fit <- jags.fit.list
load("./trunk/results/CarriageDuration/jags.fits.CASES_NOTREAT.RData")
case.notreat.jags.fit <- jags.fit.list
load("./trunk/results/CarriageDuration/jags.fits.CASES_ATB.RData")
case.atb.jags.fit <- jags.fit.list

# jags.fit.list <- list(jags.fit=jags.fit, jags.chains=jags.chains, 
#                       jags.setup=jags.setup, jags.params=jags.setup)

# print(gelman.diag(full.fit.LWW.case))
# print(gelman.diag(full.fit.LWW.carr))



# Weibull distribution parameters and key quantiles of the distribution of time to bacteria clearance for diphtheria bacillus infection.


# Carriers ...................................................
#in jags alpha = shape in R and beta = rate in R

chains_ <- carr.jags.fit$jags.chains
jags.fit_ <- c(median(chains_$v_v), quantile(chains_$v_v, prob=c(0.025,0.975)))
jags.fit_ <- rbind(jags.fit_, c(median(chains_$scale_v), quantile(chains_$scale_v, prob=c(0.025,0.975))))
for (q in c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) {
  tmp <- qweibull(q, chains_$v_v, chains_$scale_v)
  jags.fit_ <- rbind(jags.fit_, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}

tmp <- chains_$scale_v * gamma(1+1/chains_$v_v)
jags.fit_ <- rbind(quantile(tmp, prob=c(.5,.025,.975)), jags.fit_)
colnames(jags.fit_) <- c("est","CIlow","CIhigh")
rownames(jags.fit_) <- c("mean", "shape", "scale", "p5","p25","p50","p75","p95","p99")
kable(jags.fit_, format="markdown", digits=2)

carr.jags.summary <- jags.fit_


# Cases - No treatment ...................................................
#in jags alpha = shape in R and beta = rate in R

chains_ <- case.notreat.jags.fit$jags.chains
jags.fit_ <- c(median(chains_$v_v), quantile(chains_$v_v, prob=c(0.025,0.975)))
jags.fit_ <- rbind(jags.fit_, c(median(chains_$scale_v), quantile(chains_$scale_v, prob=c(0.025,0.975))))
for (q in c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) {
  tmp <- qweibull(q, chains_$v_v, chains_$scale_v)
  jags.fit_ <- rbind(jags.fit_, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}

tmp <- chains_$scale_v * gamma(1+1/chains_$v_v)
jags.fit_ <- rbind(quantile(tmp, prob=c(.5,.025,.975)), jags.fit_)
colnames(jags.fit_) <- c("est","CIlow","CIhigh")
rownames(jags.fit_) <- c("mean", "shape", "scale", "p5","p25","p50","p75","p95","p99")
kable(jags.fit_, format="markdown", digits=2)

case.notreat.jags.summary <- jags.fit_


# Cases - Antibiotics ...................................................
#in jags, alpha = shape in R and beta = rate in R

chains_ <- case.atb.jags.fit$jags.chains
jags.fit_ <- c(median(chains_$v_v), quantile(chains_$v_v, prob=c(0.025,0.975)))
jags.fit_ <- rbind(jags.fit_, c(median(chains_$scale_v), quantile(chains_$scale_v, prob=c(0.025,0.975))))
for (q in c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) {
  tmp <- qweibull(q, chains_$v_v, chains_$scale_v)
  jags.fit_ <- rbind(jags.fit_, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}

tmp <- chains_$scale_v * gamma(1+1/chains_$v_v)
jags.fit_ <- rbind(quantile(tmp, prob=c(.5,.025,.975)), jags.fit_)
colnames(jags.fit_) <- c("est","CIlow","CIhigh")
rownames(jags.fit_) <- c("mean", "shape", "scale", "p5","p25","p50","p75","p95","p99")
kable(jags.fit_, format="markdown", digits=2)

case.atb.jags.summary <- jags.fit_






# Clearance by X Day ------------------------------------------------------
days.clearance <- c(7,14,21,28,40,50,55,60,65,100)
case.notreat.dur.days <- case.atb.dur.days <- carr.dur.days <- NULL
for (d in days.clearance) {
  tmp <- pweibull(d, case.notreat.jags.fit$jags.chains$v_v, case.notreat.jags.fit$jags.chains$scale_v)
  case.notreat.dur.days <- rbind(case.notreat.dur.days, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
  
  tmp <- pweibull(d, carr.jags.fit$jags.chains$v_v, carr.jags.fit$jags.chains$scale_v)
  carr.dur.days <- rbind(carr.dur.days, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
  
  tmp <- pweibull(d, case.atb.jags.fit$jags.chains$v_v, case.atb.jags.fit$jags.chains$scale_v)
  case.atb.dur.days <- rbind(case.atb.dur.days, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}
colnames(carr.dur.days) <- colnames(case.notreat.dur.days) <- colnames(case.atb.dur.days) <- c('median', '2.5%', '97.5%')
row.names(carr.dur.days) <- row.names(case.notreat.dur.days) <- row.names(case.atb.dur.days) <- days.clearance

carr.dur.days
case.notreat.dur.days
case.atb.dur.days

1-carr.dur.days
1-case.notreat.dur.days
1-case.atb.dur.days


# Clearance percentiles ------------------------------------------------------

case.notreat.dur.days <- case.atb.dur.days <- carr.dur.days <- NULL
for (p in c(.95, .99)) {
  tmp <- qweibull(p, case.notreat.jags.fit$jags.chains$v_v, case.notreat.jags.fit$jags.chains$scale_v)
  case.notreat.dur.days <- rbind(case.notreat.dur.days, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
  
  tmp <- qweibull(p, case.atb.jags.fit$jags.chains$v_v, case.atb.jags.fit$jags.chains$scale_v)
  case.atb.dur.days <- rbind(case.atb.dur.days, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
  
  tmp <- qweibull(p, carr.jags.fit$jags.chains$v_v, carr.jags.fit$jags.chains$scale_v)
  carr.dur.days <- rbind(carr.dur.days, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}
colnames(case.notreat.dur.days) <- colnames(case.atb.dur.days) <- colnames(carr.dur.days) <- c('median', '2.5%', '97.5%')
row.names(case.notreat.dur.days) <- row.names(case.atb.dur.days) <- row.names(carr.dur.days) <- c(.95, .99)

carr.dur.days
case.notreat.dur.days
case.atb.dur.days






# RATIO of DURATIONS ------------------------------------------------------

n.iter <- 1000
carr.samples <- NULL
for (i in 1:nrow(carr.jags.fit$jags.chains)) {
  tmp <- rweibull(n.iter, carr.jags.fit$jags.chains$v_v[i], carr.jags.fit$jags.chains$scale_v[i])
  carr.samples <- c(carr.samples, as.numeric(tmp))
}

case.samples <- NULL
for (i in 1:nrow(case.notreat.jags.fit$jags.chains)) {
  tmp <- rweibull(n.iter, case.notreat.jags.fit$jags.chains$v_v[i], case.notreat.jags.fit$jags.chains$scale_v[i])
  case.samples <- c(case.samples, as.numeric(tmp))
}

case.atb.samples <- NULL
for (i in 1:nrow(case.atb.jags.fit$jags.chains)) {
  tmp <- rweibull(n.iter, case.atb.jags.fit$jags.chains$v_v[i], case.atb.jags.fit$jags.chains$scale_v[i])
  case.atb.samples <- c(case.atb.samples, as.numeric(tmp))
}


mean(carr.samples)
quantile(carr.samples, probs=c(0.025, .25, .5, .75, .975))
summary(carr.samples)

mean(case.samples)
quantile(case.samples, probs=c(0.025, .25, .5, .75, .975))
summary(case.samples)

mean(case.atb.samples)
quantile(case.atb.samples, probs=c(0.025, .25, .5, .75, .975))
summary(case.atb.samples)

case.carr.mean.ratio <- case.carr.median.ratio <- numeric(1000)
case.caseatb.mean.ratio <- case.caseatb.median.ratio <- numeric(1000)
case.caseatb.mean.diff <- case.caseatb.median.diff <- numeric(1000)

for (i in 1:1000) {
  case.tmp <- sample(case.samples, 1000, replace=TRUE)
  carr.tmp <- sample(carr.samples, 1000, replace=TRUE)
  case.atb.tmp <- sample(case.atb.samples, 1000, replace=TRUE)
  
  case.carr.mean.ratio[i] <- mean(case.tmp) / mean(carr.tmp)
  case.carr.median.ratio[i] <- median(case.tmp) / median(carr.tmp)
  
  case.caseatb.mean.ratio[i] <- mean(case.tmp) / mean(case.atb.tmp)
  case.caseatb.median.ratio[i] <- median(case.tmp) / median(case.atb.tmp)
  
  case.caseatb.mean.diff[i] <- mean(case.tmp) - mean(case.atb.tmp)
  case.caseatb.median.diff[i] <- median(case.tmp) - median(case.atb.tmp)
}
summary(case.carr.mean.ratio)
summary(case.carr.median.ratio)
quantile(case.carr.mean.ratio, probs=c(0.025, .25, .5, .75, .975))
quantile(case.carr.median.ratio, probs=c(0.025, .25, .5, .75, .975))

summary(case.caseatb.mean.ratio)
summary(case.caseatb.median.ratio)
quantile(case.caseatb.mean.ratio, probs=c(0.025, .25, .5, .75, .975))
quantile(case.caseatb.median.ratio, probs=c(0.025, .25, .5, .75, .975))

mean(case.caseatb.mean.ratio) - 1.96*sd(case.caseatb.mean.ratio)
mean(case.caseatb.mean.ratio) + 1.96*sd(case.caseatb.mean.ratio)

summary(case.caseatb.mean.diff)
quantile(case.caseatb.mean.diff, probs=c(0.025, .25, .5, .75, .975))
mean(case.caseatb.mean.diff) - 1.96*sd(case.caseatb.mean.diff)
mean(case.caseatb.mean.diff) + 1.96*sd(case.caseatb.mean.diff)




# Plot the results --------------------------------------------------------

#first make data frames to hold everything
clearance.curve.case.atb <- NULL
clearance.curve.case.notreat <- NULL
clearance.curve.carr <- NULL

for (d in seq(0,100,.1)) {
  tmp <- 1-pweibull(d, case.atb.jags.fit$jags.chains$v_v, case.atb.jags.fit$jags.chains$scale_v)
  tmp <- quantile(tmp, prob=c(0.025, .5, 0.975))
  clearance.curve.case.atb <- rbind(clearance.curve.case.atb,
                                    c(d=d, plow=tmp[1], pmid=tmp[2], phigh=tmp[3]))
  
  tmp <- 1-pweibull(d, case.notreat.jags.fit$jags.chains$v_v, case.notreat.jags.fit$jags.chains$scale_v)
  tmp <- quantile(tmp, prob=c(0.025, .5, 0.975))
  clearance.curve.case.notreat <- rbind(clearance.curve.case.notreat, 
                                        c(d=d, plow=tmp[1], pmid=tmp[2], phigh=tmp[3]))
  
  tmp <- 1-pweibull(d, carr.jags.fit$jags.chains$v_v, carr.jags.fit$jags.chains$scale_v)
  tmp <- quantile(tmp, prob=c(0.025, .5, 0.975))
  clearance.curve.carr <- rbind(clearance.curve.carr, 
                                c(d=d, plow=tmp[1], pmid=tmp[2], phigh=tmp[3]))
}

clearance.curve.case.atb <- as.data.frame(clearance.curve.case.atb)
clearance.curve.case.notreat <- as.data.frame(clearance.curve.case.notreat)
clearance.curve.carr <- as.data.frame(clearance.curve.carr)
colnames(clearance.curve.case.atb) <- colnames(clearance.curve.case.notreat) <- colnames(clearance.curve.carr) <- c("q","vcplow","vcpmid","vcphigh")
clearance.curve.case.atb$type <- 'Case, Antibiotics'
clearance.curve.case.notreat$type <- 'Case'
clearance.curve.carr$type <- 'Carrier'
clearance_curve <- as.data.frame(rbind(clearance.curve.case.atb, clearance.curve.case.notreat, clearance.curve.carr))
clearance_curve$type <- factor(clearance_curve$type, levels=c('Carrier', 'Case', 'Case, Antibiotics'))

# Save Curve Results    
save(clearance_curve, file = "./trunk/Results/carriage_curve.RData")
load(file = "./trunk/Results/carriage_curve.RData") # clearance_curve

bc.plt <- ggplot(clearance_curve, aes(x=q), group=type) +
  geom_ribbon(aes(ymin=vcplow, ymax=vcphigh, fill=type),  alpha=.4) +
  geom_line(aes(y=vcpmid, colour=type)) + theme_bw() +
  scale_fill_manual(values=(c("darkorange", "red", "gold"))) +
  scale_colour_manual(values=(c("darkorange", "red", "gold"))) +
  scale_x_continuous(limits=c(0, 70), expand = c(0, 0)) + 
  ylab("% w/ detectable bacilli") + xlab("days")
bc.plt








