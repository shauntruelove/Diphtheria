
# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('coda')) install.packages('coda'); library(coda)
if(!require('rjags')) install.packages('rjags'); library(rjags)


recalc <- F

# Check if Results Exist -------------------------------------------

if (!file.exists("./results/full.fit.jags.LWW.TREATED.RData") | recalc) {
  
  # Set up JAGS global parameters
  n.adapt <- 10000
  iters <- 1000000
  thin <- 500
  n.chains <- 3
  
  # FUNCTIONS FOR JAGS FITTING ----------------------------------------------
  
  # Function for setting up JAGS
  # Function for setting up JAGS
  carriage.jags.setup.fn <- function(data) {
    
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
    
    # add to right bound to fix problem with having left and right bounds the same
    jags.data$CR <- jags.data$CR + 0.1
    
    #set the parameters we want to track
    #  - v_v     : weibull shape
    #  - scale_v : weibull scale
    parameters <- c("v_v","scale_v", "loglik")  
    
    set.seed(12345) #if this is not included, multiple initializations maybe needed.
    #initialization function for jags
    jags.inits <-  function() {
      rc <-list(
        #E = rep(0.0000000011, n.subjects), #start with a fixed E to avoid bad starting points
        v_v = runif(1,1,10),
        scale_v = runif(1,0,1), 
        Y_V=CDyInit)
      print(rc)
      return(rc)
    }
    
    return(list(jags.data=jags.data, jags.inits=jags.inits, parameters=parameters))
  }
  
  
  # Function to fit JAGS model, clean, and save
  carriage.jags.fit.fn <- function(jagsfit.LWW=jagsfit.LWW, jags.setup=carr.jags.setup, group='CARRIERS',
                                   iters, thin, n.chains, save.jags.fit=TRUE) {
    
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
    if(save.jags.fit) {
      save(jags.fit.list, file=paste0("./results/jags.fits.",group,".RData"))
    }
    return(jags.fit.list)
  }
  

  
  
  
  # Load the raw colonization data ------------------------------------------

  dc.data <- read.csv(file = "./data/duration_colonization_treated.csv", header=TRUE)
  dc.data$Treatment <- tolower(dc.data$Treatment)
  
  # Convert to long
  dc.long <- data.frame(Case=rep(dc.data$Case, times=dc.data$Num.convert),
                        Treatment=rep(dc.data$Treatment, times=dc.data$Num.convert),
                        Antitoxin=rep(dc.data$Antitoxin, times=dc.data$Num.convert),
                        carr.min=rep(dc.data$carr.min, times=dc.data$Num.convert), 
                        carr.max=rep(dc.data$carr.max, times=dc.data$Num.convert),
                        study=rep(dc.data$Study, times=dc.data$Num.convert))
  
  
  
  
  
  #///////////////////////////////////////////////////////////////////////////////////////
  # Bayesian MCMC -----------------------------------------------------------
  #///////////////////////////////////////////////////////////////////////////////////////
  
  
  #///////////////////////////////////////////////////////////////////////////////////////
  # Shared colonization duration ----------------------------------------------------------------
  
  # Setup data and initials
  jags.setup <- carriage.jags.setup.fn(dc.long)
  
  #initialize JAGS model
  jags.model.file <- ifelse(length(jags.setup$jags.data$CDrc)==0, 
                            './source/DistributionFitCarriage_NoRC.jags', './source/DistributionFitCarriage.jags')
  jagsfit.LWW <- jags.model(file=jags.model.file, 
                            data=jags.setup$jags.data, 
                            inits=jags.setup$jags.inits, 
                            n.chains=n.chains, quiet=F, 
                            n.adapt=n.adapt)
  
  # Fit JAGS model, clean, and save
  jags.fit.list <- carriage.jags.fit.fn(jagsfit.LWW = jagsfit.LWW, jags.setup=jags.setup, 
                                        group='TREATED', iters, thin, n.chains, save.jags.fit = T)
  
}
  
  
  
# JAGS Results ------------------------------------------------------------

load("./results/jags.fits.TREATED.RData")
treated.jags.fit <- jags.fit.list


chains_ <- treated.jags.fit$jags.chains
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




# Clearance by X Day ------------------------------------------------------
days.clearance <- c(7,14,21,28,40,50,55,60,65,100)
treated.dur.days <- NULL
for (d in days.clearance) {
  tmp <- pweibull(d, treated.jags.fit$jags.chains$v_v, treated.jags.fit$jags.chains$scale_v)
  treated.dur.days <- rbind(treated.dur.days, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}
colnames(treated.dur.days) <- c('median', '2.5%', '97.5%')
row.names(treated.dur.days) <- days.clearance

treated.dur.days
1-treated.dur.days


# Clearance percentiles ------------------------------------------------------

treated.dur.days <- NULL
for (p in c(.95, .99)) {
  tmp <- qweibull(p, treated.jags.fit$jags.chains$v_v, treated.jags.fit$jags.chains$scale_v)
  treated.dur.days <- rbind(treated.dur.days, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}
colnames(treated.dur.days) <- c('median', '2.5%', '97.5%')
row.names(treated.dur.days) <- c(.95, .99)

treated.dur.days




# Plot the results --------------------------------------------------------

#first make data frames to hold everything
clearance.curve.treated <- NULL

for (d in seq(0,100,.1)) {
  tmp <- 1-pweibull(d, treated.jags.fit$jags.chains$v_v, treated.jags.fit$jags.chains$scale_v)
  tmp <- quantile(tmp, prob=c(0.025, .5, 0.975))
  clearance.curve.treated <- rbind(clearance.curve.treated, c(d=d, plow=tmp[1], pmid=tmp[2], phigh=tmp[3]))
}
clearance_curve_treated <- as.data.frame(clearance.curve.treated)


colnames(clearance_curve_treated) <- c("q","vcplow","vcpmid","vcphigh")
clearance_curve_treated$type <- 'treated'

# Save Curve Results    
# save(clearance_curve_treated, file = "results/carriage_curve_treated.RData")
# load(file = "results/carriage_curve_treated.RData") # clearance_curve_treated
# write.csv(clearance_curve_treated, file='results/clearance_curve_treated.csv', row.names = FALSE)
clearance_curve_treated <- read.csv(file='results/clearance_curve_treated.csv', header=TRUE)


bc.plt <- ggplot(clearance_curve_treated, aes(x=q), group=type) +
  geom_ribbon(aes(ymin=vcplow, ymax=vcphigh, fill=type),  alpha=.4) +
  geom_line(aes(y=vcpmid, colour=type)) + theme_bw() +
  scale_fill_manual(values=(c("darkorange"))) +
  scale_colour_manual(values=(c("darkorange"))) +
  scale_x_continuous(limits=c(0, 70), expand = c(0, 0)) + 
  ylab("% w/ detectable bacilli") + xlab("days")
bc.plt
  


# Plot comparison with no treatment

clearance_curve_shared <- read.csv(file='results/clearance_curve_shared.csv', header=TRUE)
clearance_curve_shared <- clearance_curve_shared %>% mutate(Antibiotics="No")
clearance_curve_treated <- clearance_curve_treated %>% mutate(Antibiotics="Yes")
clearance_curve_data <- rbind(clearance_curve_treated, clearance_curve_shared)


compare.plt <- ggplot(clearance_curve_data, aes(x=q, group=Antibiotics)) +
  geom_ribbon(aes(ymin=vcplow, ymax=vcphigh, fill=Antibiotics),  alpha=.4) +
  geom_line(aes(y=vcpmid, colour=Antibiotics)) + theme_bw() +
  scale_fill_manual(values=(c("darkorange", "blue"))) +
  scale_colour_manual(values=(c("darkorange", "blue"))) +
  scale_x_continuous(limits=c(0, 70), expand = c(0, 0)) + 
  ylab("% w/ detectable bacilli") + xlab("days")
compare.plt

