# Install packages
if(!require('fitdistrplus')) install.packages('fitdistrplus'); library(fitdistrplus)
if(!require('coda')) install.packages('coda'); library(coda)
if(!require('rjags')) install.packages('rjags'); library(rjags)
if(!require('lubridate')) install.packages('lubridate'); library(lubridate)

if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)
if(!require('knitr')) install.packages('knitr'); library(knitr)

recalc <- F

# Check if Results Exist -------------------------------------------

if (!file.exists("results/full.fit.jags.LWW.Incub.RData") | recalc) {
    
  
    incub.data <- read.csv(file='./data/incubation_period.csv', header=TRUE)
    
    # Normalize the times and give range within the day
    incub.data$ER[is.na(incub.data$ER) & !is.na(incub.data$SL)] <- incub.data$SL[is.na(incub.data$ER) & !is.na(incub.data$SL)]
    incub.data[,c('EL','ER','SL','SR')] <- apply(incub.data[,c('EL','ER','SL','SR')], 2, lubridate::mdy)
    incub.data$ER <- incub.data$ER+(23/24)
    incub.data$SR <- incub.data$SR+(23/24)
    incub.data$ER[incub.data$ER>=incub.data$SL] <- incub.data$SL[incub.data$ER>=incub.data$SL] - 1/24
    
    unique(incub.data$Study)
    
    # Recalculate incubation period to days since S1
    time0 <- incub.data[,'EL']
    incub.data[,c('EL','ER','SL','SR')] <- incub.data[,c('EL','ER','SL','SR')] - time0
    
    #incub.data <- read.csv('./data/zika.anal.jags.csv', header=T)
    
    #make a data object for JAGS
    jags.data <- list(
      #EL=incub.data.jags$EL,
      ER=incub.data$ER,
      SL=incub.data$SL,
      SR=incub.data$SR)
    
    #Let jags know the censoring situation
    # 1 means interval censored
    # 2 means event occured after known time
    jags.data$IPisCensored <- rep(1, length(jags.data$ER))
    jags.data$IPisCensored[which(is.na(jags.data$SR))] <- 2

    # number of subjects
    n.subjects <- nrow(incub.data)
    
    #define variable to hold the length of the time to event for diphtheria clearance
    jags.data$Y_S <- rep(NA, n.subjects)
    
    #set the intitial values for time to event (i.e., Y_S)
    IPyInit <- as.numeric(jags.data$SL)
    IPyInit[which(is.na(jags.data$SL)==T)] <- 0
    IPyInit[which(IPyInit==0)] <- 0.0000000001
    
    #identify which observations are interval censored and which are right censored. 
    jags.data$dic <- which(jags.data$IPisCensored==1) #this is all inclubation period observtions
    
    #set the parameters we want to track
    parameters <- c("lm","lsd","E")
    
    set.seed(12345) #if this is not included, multiple initializations maybe needed.
    #initialization function for jags
    jags.inits <-  function() {
        rc <-list(E=rep(0.0000000000011,23), #start with a fixed E to avoid bad starting points
                  lm=runif(1,log(2),log(10)),
                  lsd=runif(1,.1,log(3)),
                  Y_S=IPyInit)
      print(rc)
      return(rc)
    }
    
    #initialize JAGS model
    jagsfit.LWW <- jags.model(file='./source/DistributionFitIncub_lognormal.jags', 
                              data=jags.data, 
                              inits=jags.inits, 
                              n.chains=3, quiet=F, 
                              n.adapt=10000)
    iters <- 1000000
    thin  <- 500
    
    full.fit.LWW <- coda.samples(jagsfit.LWW, parameters, n.iter=iters, thin=thin, n.chains=3)
    
    #make all of the chains a single matrix with a burnin removed
    ABC1=as.matrix(full.fit.LWW[[1]][,])
    ABC2=as.matrix(full.fit.LWW[[2]][,])
    ABC3=as.matrix(full.fit.LWW[[3]][,])
    
    ABC1=ABC1[(thin+1):(iters/thin),]
    ABC2=ABC2[(thin+1):(iters/thin),]
    ABC3=ABC3[(thin+1):(iters/thin),]

    #recreate MCMC object for diagnostics
    full.fit.LWW.i<-list(as.mcmc(ABC1), as.mcmc(ABC2), as.mcmc(ABC3)) 
    chains.LWW.i <- rbind(ABC1,ABC2,ABC3)
    colnames(chains.LWW.i) <- varnames(full.fit.LWW.i[[1]])
    chains.LWW.i <- as.data.frame(chains.LWW.i)
    save(full.fit.LWW.i, chains.LWW.i, file="./results/full.fit.jags.LWW.Incub.LogNormal.RData")
}




# JAGS Results ------------------------------------------------------------

load("./results/full.fit.jags.LWW.Incub.LogNormal.RData")  # loads full.fit.LWW.i and chains.LWW.i



# Incubation Period ---------------------------------------------------------

#in stan alpha = shape in R and beta = rate in R
bshed.fit.jags.LWW.i <- c(median(chains.LWW.i$lm), quantile(chains.LWW.i$lm, prob=c(0.025,0.975)))
bshed.fit.jags.LWW.i <- rbind(bshed.fit.jags.LWW.i, c(median(chains.LWW.i$lsd), quantile(chains.LWW.i$lsd, prob=c(0.025,0.975))))
for (q in c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99)) {
    tmp <- qlnorm(q, chains.LWW.i$lm, chains.LWW.i$lsd)
    bshed.fit.jags.LWW.i <- rbind(bshed.fit.jags.LWW.i, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}

tmp <- exp(chains.LWW.i$lm + ((chains.LWW.i$lsd^2)/2))
bshed.fit.jags.LWW.i <- rbind(quantile(tmp, prob=c(.5,.025,.975)), bshed.fit.jags.LWW.i)
colnames(bshed.fit.jags.LWW.i) <- c("est","CIlow","CIhigh")
rownames(bshed.fit.jags.LWW.i) <- c("mean",
                                  "meanlog",
                                  "sdlog",
                                  "p2.5", "p5","p25","p50","p75","p95","p97.5","p99")
kable(bshed.fit.jags.LWW.i, format="markdown", digits=2)




# Plot the results --------------------------------------------------------

# data frames to hold everything
Incub.curve.i.lnorm <- NULL
options(scipen=999)

for (d in seq(0,20,.1)) {
    tmp.lnorm <- dlnorm(d, meanlog=chains.LWW.i$lm, sdlog=chains.LWW.i$lsd)
    tmp.lnorm <- quantile(tmp.lnorm, prob=c(0.025, .5, 0.975))
    Incub.curve.i.lnorm <- rbind(Incub.curve.i.lnorm, c(d=d,
                                                      plow=tmp.lnorm[1],
                                                      pmid=tmp.lnorm[2],
                                                      phigh=tmp.lnorm[3]))
}

Incub.curve.i <- as.data.frame(Incub.curve.i.lnorm)
Incub.curve.i$type <- 'Incubation'
colnames(Incub.curve.i) <- c("q","plow","pmid","phigh","type")


save(Incub.curve.i, file = "results/Incub_curve_lnorm.Rdata")
# Save Curve data as CSV
#write.csv(Incub.curve.i, file='results/incubation_curve.csv', row.names = FALSE)

load(file = "results/Incub_curve_lnorm.Rdata") # Incub.curve.i
# Read CSV
Incub.curve.i <- read.csv(file='results/incubation_curve.csv', header=TRUE)


# Plot Incubation Period

incub.plot <- ggplot(Incub.curve.i, aes(x=q)) +
    geom_ribbon(aes(ymin=plow, ymax=phigh, fill=type),  alpha=.4) +
    geom_line(aes(y=pmid, colour=type, group=type)) + theme_bw() +
    scale_fill_manual(values="green") +
    scale_colour_manual(values="green") +
    scale_x_continuous(limits=c(0, 20), expand = c(0, 0)) + 
    ylab("p") + xlab("days")
incub.plot

