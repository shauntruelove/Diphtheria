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

if (!file.exists("./trunk/results/full.fit.jags.LWW.GenTime.RData") | recalc) {

    # Load the data
    gentime.data <- read.csv(file='./data/generation_time.csv', header = TRUE, stringsAsFactors = F)
    
    
    gentime.data[,c('S1L','S1R','S2L','S2R')] <- apply(gentime.data[,c('S1L','S1R','S2L','S2R')], 2, lubridate::mdy)
    #gentime.dat[,c('S1L','S1R','S2L','S2R')] <- apply(gentime.dat[,c('S1L','S1R','S2L','S2R')], 2, as.Date, "%m-%d-%Y")
    unique(gentime.data$Study)
    
    
    # Recalculate generation time to days since S1
    time0 <- gentime.data[,'S1L']
    gentime.data[,c('S1L','S1R','S2L','S2R')] <- gentime.data[,c('S1L','S1R','S2L','S2R')] - time0

    # Adjust the minimum times
    gentime.data <- gentime.data %>%
        mutate(GTL = S2L-S1R,
               GTR = S2R) %>% as.data.frame
    gentime.data$GTL <- gentime.data$GTL - 23/24   # Could be 1-23 hours during the day
    gentime.data$GTL[gentime.data$GTL<=0] <- 1/24  # Set min to 1 hr
    gentime.data$GTR <- gentime.data$GTR + 23/24   # Could be 23rd hour of the day
    
    
    # Bayesian MCMC -----------------------------------------------------------
    
    #assuming observations with no upper limite have NAs in the R slot. Currently only in carr.max column
    gentime.data.jags <- gentime.data[,c('ID', 'GTL', 'GTR')]
    gentime.data.jags$GTR[which(is.infinite(gentime.data.jags$GTR))] <-  NA
    
    #missing observation are just unbounded, replace with 0
    gentime.data.jags$GTL[which(is.na(gentime.data.jags$GTL))] <-  0
    
    
    #make a data object for JAGS
    jags.data <- list(
        GTL=gentime.data.jags$GTL,
        GTR=gentime.data.jags$GTR)
    
    
    #Let jags know the censoring situation
    # 1 means interval censored
    # 2 means event occured after known time
    jags.data$VSPisCensored=rep(1, length(jags.data$GTL))
    jags.data$VSPisCensored[which(is.na(jags.data$GTR))]=2
    
    
    n.subjects <- nrow(gentime.data.jags)
    
    #define variable to hold the length of the time to event for bacilli clearance
    jags.data$Y_GT <- rep(NA, n.subjects)
    
    #set the intitial values for time to event (i.e., Y_GT)
    VSPyInit <- jags.data$GTL
    VSPyInit[which(is.na(jags.data$GTL)==T)] <- 0
    VSPyInit[which(VSPyInit==0)] <- 0.000000001
    
    #identify which observations are interval censored and which are right censored. 
    jags.data$GTdic <- which(jags.data$VSPisCensored==1)
    #jags.data$GTrc  <- which(jags.data$VSPisCensored==2)
    
    # # Make GTrc NA if empty
    # if (length(jags.data$GTrc)==0) jags.data$GTrc <- NA
    # 
    #set the parameters we want to track
    parameters <- c("shape_gt","scale_gt")
    
    set.seed(12345) #if this is not included, multiple initializations maybe needed.
    #initialization function for jags
    jags.inits <-  function() {
        rc <-list(
                  shape_gt = runif(1,0,20),
                  scale_gt = runif(1,0,10), 
                  Y_GT=VSPyInit)
        print(rc)
        return(rc)
    }
    
    
    
    #not evaluated during knitting, run to get the next section to work
    #note, we deal with the burnin manually later.
    
    #initialize JAGS model
    jagsfit.LWW <- jags.model(file='./source/DistributionFitGenTime_gamma.jags', 
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
    # ABC1=ABC1[5001:(iters/thin),]
    # ABC2=ABC2[5001:(iters/thin),]
    # ABC3=ABC3[5001:(iters/thin),]
    
    #recreate MCMC object for diagnostics
    full.fit.LWW.gt<-list(as.mcmc(ABC1), as.mcmc(ABC2), as.mcmc(ABC3)) 
    chains.LWW.gt <- rbind(ABC1,ABC2,ABC3)
    colnames(chains.LWW.gt) <- varnames(full.fit.LWW.gt[[1]])
    chains.LWW.gt <- as.data.frame(chains.LWW.gt)
    save(full.fit.LWW.gt, chains.LWW.gt, file="./results/full.fit.jags.LWW.GenTime.RData")
    
}




# JAGS Results ------------------------------------------------------------

load("./results/full.fit.jags.LWW.GenTime.RData")



# Table S3: Weibull distribution parameters and key quantiles of the distribution of time to bacteria clearance for diphtheria bacillus infection.

#in stan alpha = shape in R and beta = rate in R
bshed.fit.jags.LWW.gt <- c(median(chains.LWW.gt$shape_gt), quantile(chains.LWW.gt$shape_gt,prob=c(0.025,0.975)))
bshed.fit.jags.LWW.gt <- rbind(bshed.fit.jags.LWW.gt, c(median(chains.LWW.gt$scale_gt), quantile(chains.LWW.gt$scale_gt,prob=c(0.025,0.975))))
for (q in c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) {
    tmp <- qgamma(q, chains.LWW.gt$shape_gt, scale=chains.LWW.gt$scale_gt)
    bshed.fit.jags.LWW.gt <- rbind(bshed.fit.jags.LWW.gt, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}

tmp <- chains.LWW.gt$scale_gt * chains.LWW.gt$shape_gt
bshed.fit.jags.LWW.gt <- rbind(quantile(tmp, prob=c(.5,.025,.975)), bshed.fit.jags.LWW.gt)
colnames(bshed.fit.jags.LWW.gt) <- c("est","CIlow","CIhigh")
rownames(bshed.fit.jags.LWW.gt) <- c("mean",
                                  "shape",
                                  "scale",
                                  "p5","p25","p50","p75","p95","p99")
kable(bshed.fit.jags.LWW.gt, format="markdown", digits=2)

# variance
scale_gt <- bshed.fit.jags.LWW.gt[3,1]
shape_gt <- bshed.fit.jags.LWW.gt[2,1]
var_gt <- scale_gt^2 * shape_gt
var_gt


# Generation Time Parameters

gentime_params <- list(shape_mean=mean(chains.LWW.gt$shape_gt), shape_sd=sd(chains.LWW.gt$shape_gt), 
                       scale_mean=mean(chains.LWW.gt$scale_gt), scale_sd=sd(chains.LWW.gt$scale_gt))
save(gentime_params, file="./results/gentime_params.RData")





# Gneration time by X Day ------------------------------------------------------

GT.dist <- NULL
for (t in c(.25,.5,.75,1:10, 15,20)) {
    tmp <- pgamma(t, chains.LWW.gt$shape_gt, scale=chains.LWW.gt$scale_gt) 
    GT.dist <- rbind(GT.dist, c(median(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}
colnames(GT.dist) <- c('median', '2.5%', '97.5%')
row.names(GT.dist) <- c(.25,.5,.75,1:10, 15,20)


GT.dist
1-GT.dist




# Plot the results --------------------------------------------------------

# data frames to hold everything
GT.curve.gt <- NULL
options(scipen=999)

for (d in seq(0,20,.1)) {

    tmp <- dgamma(d, shape=chains.LWW.gt$shape_gt, scale=chains.LWW.gt$scale_gt)
    tmp <- quantile(tmp, prob=c(0.025, .5, 0.975))
    GT.curve.gt <- rbind(GT.curve.gt, c(d=d,
                                  plow=tmp[1],
                                  pmid=tmp[2],
                                  phigh=tmp[3]))
    
}


GT.curve.gt <- as.data.frame(GT.curve.gt)
colnames(GT.curve.gt) <- c("q","plow","pmid","phigh")
GT.curve.gt$type <- 'Generation Time'
GT_curve <- GT.curve.gt
save(GT_curve, file = "results/GenTime_curve.Rdata")

# Save Curve data as CSV
write.csv(GT_curve, file='results/serial_interval_curve.csv', row.names = FALSE)


load(file = "results/GenTime_curve.Rdata") # GT_curve
# Read CSV
GT_curve <- read.csv(file='results/serial_interval_curve.csv', header=TRUE)


# Plot Generation Time
gt.plot <- ggplot(GT_curve, aes(x=q)) +
    geom_ribbon(aes(ymin=plow, ymax=phigh, fill=type),  alpha=.4) +
    geom_line(aes(y=pmid, colour=type)) + theme_bw() +
    scale_fill_manual(values="blue") +
    scale_colour_manual(values="blue") +
    scale_x_continuous(limits=c(0, 20), expand = c(0, 0)) + 
    ylab("p") + xlab("days")
gt.plot

