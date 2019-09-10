

# Load Required Packages --------------------------------------------------

if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('rstan')) install.packages('rstan'); library(rstan)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('scales')) install.packages('scales'); library(scales)
if(!require('reshape2')) install.packages('reshape2'); library(scales)

 

# Source files ------------------------------------------------------------

source('source/Vc_and_R0_functions.R')



data <- vc_treat(delay=1, alpha=0, restrict_quantiles=c(0,1))
R <- numeric(nrow(data))

for (i in 1:nrow(data)) {
  R[i] <- R_calc(V3=0, Rs=data$Rs[i], p0=data$p0[i], p3=data$p3[i], tau=data$tau[i], delta=calc_delta(delay=1), alpha=0.90)
}

summary(R)
quantile(R, probs=c(0.025,0.975))





# Setup Stan Data ---------------------------------------------------------

# Get the stan results and restrict to the 95% CI of R0
R095 <- quantile(fits$R0_preds, probs=c(0.025,0.975))
fits.df <- data.frame(R0_preds=fits$R0_preds, lRs_star=fits$lRs_star, Rs=fits$Rs_tmp, pk=fits$pk, tau=fits$tau, Vc=fits$Vc)
fits.df$R0 <- fits.df$pk.1*fits.df$Rs + (1-fits.df$pk.1)*fits.df$tau*fits.df$Rs
fits.df <- fits.df %>% rename(p0=pk.1, p2=pk.2, p3=pk.3) %>% as.data.frame()
fits95 <- fits.df %>% filter(R0>=R095[1] & R0<=R095[2])

  


# Calculate treatment-adjusted R for each iteration -----------------------

V3 <- seq(0,1,.01)

R_0pct <- matrix(NA, nrow=nrow(fits95), ncol=length(V3))
colnames(R_0pct) <- V3
row.names(R_0pct) <- 1:nrow(fits95)
R_25pct <- R_50pct <- R_75pct <- R_100pct <- R_0pct
for (i in 1:nrow(fits95)) {
    R_0pct[i,] <- R_calc(V3=V3, Rs=fits95$Rs[i], p0=fits95$p0[i], p3=fits95$p3[i], tau=fits95$tau[i], delta=delta.2day, alpha=0)
    R_25pct[i,] <- R_calc(V3=V3, Rs=fits95$Rs[i], p0=fits95$p0[i], p3=fits95$p3[i], tau=fits95$tau[i], delta=delta.2day, alpha=.25)
    R_50pct[i,] <- R_calc(V3=V3, Rs=fits95$Rs[i], p0=fits95$p0[i], p3=fits95$p3[i], tau=fits95$tau[i], delta=delta.2day, alpha=.50)
    R_75pct[i,] <- R_calc(V3=V3, Rs=fits95$Rs[i], p0=fits95$p0[i], p3=fits95$p3[i], tau=fits95$tau[i], delta=delta.2day, alpha=.75)
    R_100pct[i,] <- R_calc(V3=V3, Rs=fits95$Rs[i], p0=fits95$p0[i], p3=fits95$p3[i], tau=fits95$tau[i], delta=delta.2day, alpha=1.0)
}

R_0tmp <- data.frame(V3, t(R_0pct))
R_0pct_long <- melt(R_0tmp, id.var="V3")
R_25tmp <- data.frame(V3, t(R_25pct))
R_25pct_long <- melt(R_25tmp, id.var="V3")
R_50tmp <- data.frame(V3, t(R_50pct))
R_50pct_long <- melt(R_50tmp, id.var="V3")
R_75tmp <- data.frame(V3, t(R_75pct))
R_75pct_long <- melt(R_75tmp, id.var="V3")
R_100tmp <- data.frame(V3, t(R_100pct))
R_100pct_long <- melt(R_100tmp, id.var="V3")

R_all_treat <- rbind(data.frame(R_0pct_long, alpha=0), data.frame(R_25pct_long, alpha=.25), data.frame(R_50pct_long, alpha=.50),
                     data.frame(R_75pct_long, alpha=.75), data.frame(R_100pct_long, alpha=1))
R_all_treat <- rbind(data.frame(R_0pct_long, alpha=0), data.frame(R_50pct_long, alpha=.50),
                     data.frame(R_100pct_long, alpha=1))
colnames(R_all_treat) <- c("V3", "Stan_iter", "R", "alpha")


# Get quantiles of the results too
R_0_quant <- R_25_quant <- R_50_quant <- R_75_quant <- R_100_quant <- data.frame(V3, LowCI=NA, q25=NA, Median=NA, q75=NA, HighCI=NA, Mean=NA)
for (i in 1:length(V3)){
    R_0_quant[i,-1] <- c(quantile(R_0pct[,i], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(R_0pct[,i]))
    R_25_quant[i,-1] <- c(quantile(R_25pct[,i], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(R_25pct[,i]))
    R_50_quant[i,-1] <- c(quantile(R_50pct[,i], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(R_50pct[,i]))
    R_75_quant[i,-1] <- c(quantile(R_75pct[,i], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(R_75pct[,i]))
    R_100_quant[i,-1] <- c(quantile(R_100pct[,i], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(R_100pct[,i]))
}

R_quants <- rbind(data.frame(alpha=0, R_0_quant), data.frame(alpha=.25, R_25_quant), data.frame(alpha=.50, R_50_quant),
                  data.frame(alpha=.75,R_75_quant), data.frame(alpha=100, R_100_quant))
R_quants_long <- melt(R_quants, id.vars = c("V3", "alpha"))
colnames(R_quants_long) <- c("V3", "alpha", "Variable", "R")

R_quants_2daydelay <- R_quants
R_quants_long_2daydelay <- R_quants_long

save(R_quants, R_quants_long, R_all_treat, file="./trunk/results/R_by_V3-2daydelay.RData")




# Calculate treatment-adjusted R for each iteration - 0 day delay -----------------------

V3 <- seq(0,1,.01)

R_0pct <- matrix(NA, nrow=nrow(fits95), ncol=length(V3))
colnames(R_0pct) <- V3
row.names(R_0pct) <- 1:nrow(fits95)
R_25pct <- R_50pct <- R_75pct <- R_100pct <- R_0pct
for (i in 1:nrow(fits95)) {
    R_0pct[i,] <- R_calc(V3=V3, Rs=fits95$Rs[i], p0=fits95$p0[i], p3=fits95$p3[i], tau=fits95$tau[i], delta=delta.0day, alpha=0)
    R_25pct[i,] <- R_calc(V3=V3, Rs=fits95$Rs[i], p0=fits95$p0[i], p3=fits95$p3[i], tau=fits95$tau[i], delta=delta.0day, alpha=.25)
    R_50pct[i,] <- R_calc(V3=V3, Rs=fits95$Rs[i], p0=fits95$p0[i], p3=fits95$p3[i], tau=fits95$tau[i], delta=delta.0day, alpha=.50)
    R_75pct[i,] <- R_calc(V3=V3, Rs=fits95$Rs[i], p0=fits95$p0[i], p3=fits95$p3[i], tau=fits95$tau[i], delta=delta.0day, alpha=.75)
    R_100pct[i,] <- R_calc(V3=V3, Rs=fits95$Rs[i], p0=fits95$p0[i], p3=fits95$p3[i], tau=fits95$tau[i], delta=delta.0day, alpha=1.0)
}

R_0tmp <- data.frame(V3, t(R_0pct))
R_0pct_long <- melt(R_0tmp, id.var="V3")
R_25tmp <- data.frame(V3, t(R_25pct))
R_25pct_long <- melt(R_25tmp, id.var="V3")
R_50tmp <- data.frame(V3, t(R_50pct))
R_50pct_long <- melt(R_50tmp, id.var="V3")
R_75tmp <- data.frame(V3, t(R_75pct))
R_75pct_long <- melt(R_75tmp, id.var="V3")
R_100tmp <- data.frame(V3, t(R_100pct))
R_100pct_long <- melt(R_100tmp, id.var="V3")

R_all_treat <- rbind(data.frame(R_0pct_long, alpha=0), data.frame(R_25pct_long, alpha=.25), data.frame(R_50pct_long, alpha=.50),
                     data.frame(R_75pct_long, alpha=.75), data.frame(R_100pct_long, alpha=1))
R_all_treat <- rbind(data.frame(R_0pct_long, alpha=0), data.frame(R_50pct_long, alpha=.50),
                     data.frame(R_100pct_long, alpha=1))
colnames(R_all_treat) <- c("V3", "Stan_iter", "R", "alpha")


# Get quantiles of the results too
R_0_quant <- R_25_quant <- R_50_quant <- R_75_quant <- R_100_quant <- data.frame(V3, LowCI=NA, q25=NA, Median=NA, q75=NA, HighCI=NA, Mean=NA)
for (i in 1:length(V3)){
    R_0_quant[i,-1] <- c(quantile(R_0pct[,i], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(R_0pct[,i]))
    R_25_quant[i,-1] <- c(quantile(R_25pct[,i], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(R_25pct[,i]))
    R_50_quant[i,-1] <- c(quantile(R_50pct[,i], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(R_50pct[,i]))
    R_75_quant[i,-1] <- c(quantile(R_75pct[,i], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(R_75pct[,i]))
    R_100_quant[i,-1] <- c(quantile(R_100pct[,i], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)), mean(R_100pct[,i]))
}

R_quants <- rbind(data.frame(alpha=0, R_0_quant), data.frame(alpha=.25, R_25_quant), data.frame(alpha=.50, R_50_quant),
                  data.frame(alpha=.75,R_75_quant), data.frame(alpha=100, R_100_quant))
R_quants_long <- melt(R_quants, id.vars = c("V3", "alpha"))
colnames(R_quants_long) <- c("V3", "alpha", "Variable", "R")

R_quants_0daydelay <- R_quants
R_quants_long_0daydelay <- R_quants_long

save(R_quants, R_quants_long, R_all_treat, file="./trunk/results/R_by_V3-0daydelay.RData")




# Plot Scatter Plot -------------------------------------------------------

R_all_treat_samp <- R_all_treat[sample(1:nrow(R_all_treat), 5000),]
R_all_treat_samp$alpha <- factor(R_all_treat_samp$alpha)
p <- ggplot(R_all_treat_samp, aes(x=V3, y=R, group=alpha, color=alpha)) + 
    geom_point(aes(colour = alpha)) + stat_smooth(aes(colour = alpha)) +
    theme_bw()
p



# Plot CIs ----------------------------------------------------------------

R_quants <- R_quants %>% filter(alpha != 0.25, alpha !=.75)
R_quants$alpha <- factor(R_quants$alpha)


R_quants_long<- R_quants_long %>% filter(alpha != 0.25, alpha !=.75)
R_quants_long$alpha <- factor(R_quants_long$alpha)

R_quants_med <- R_quants_long %>% filter(Variable=='Median')
R_quants_CIs <- R_quants_long %>% filter(Variable %in% c("LowCI", "HighCI"))


p2 <- ggplot(R_quants, aes(x=V3, y=Median, color=alpha, group=alpha)) + geom_smooth(aes(color=alpha)) + 
    geom_hline(yintercept = 1, color='red', size=1.5) +
    geom_ribbon(aes(x=V3, ymin=LowCI, ymax=HighCI, fill=alpha, group=alpha), alpha=.4) +
    xlab("V3") + ylab("R") +
    theme_bw()
p2


p3 <- ggplot(R_quants, aes(x=V3, y=Median, color=alpha, group=alpha)) + geom_smooth(aes(color=alpha)) + 
    geom_hline(yintercept = 1, color='red', size=1.5) +
    geom_ribbon(aes(x=V3, ymin=q25, ymax=q75, fill=alpha, group=alpha), alpha=.4) +
    xlab("V3") + ylab("R") +
    theme_bw()
p3

#




# Calculate treatment-adjusted R for rohingya -----------------------------

V3 = c(0, .25, .50, .75, 1)

R_roh_0pct <- matrix(NA, nrow=nrow(fits.df), ncol=length(V3))
colnames(R_roh_0pct) <- V3
row.names(R_roh_0pct) <- 1:nrow(fits.df)
R_roh_25pct <- R_roh_50pct <- R_roh_75pct <- R_roh_100pct <- R_roh_0pct

for (i in 1:nrow(fits.df)) {
    R_roh_0pct[i,]   <- R_calc(V3=V3, Rs=fits$Rs[i,1], p0=fits$pk[i,1], p3=fits$pk[i,3], tau=fits$tau[i], delta=delta.2day, alpha=0)
    R_roh_25pct[i,]  <- R_calc(V3=V3, Rs=fits$Rs[i,1], p0=fits$pk[i,1], p3=fits$pk[i,3], tau=fits$tau[i], delta=delta.2day, alpha=.25)
    R_roh_50pct[i,]  <- R_calc(V3=V3, Rs=fits$Rs[i,1], p0=fits$pk[i,1], p3=fits$pk[i,3], tau=fits$tau[i], delta=delta.2day, alpha=.50)
    R_roh_75pct[i,]  <- R_calc(V3=V3, Rs=fits$Rs[i,1], p0=fits$pk[i,1], p3=fits$pk[i,3], tau=fits$tau[i], delta=delta.2day, alpha=.75)
    R_roh_100pct[i,] <- R_calc(V3=V3, Rs=fits$Rs[i,1], p0=fits$pk[i,1], p3=fits$pk[i,3], tau=fits$tau[i], delta=delta.2day, alpha=1.0)
}

R_roh_all <- rbind(data.frame(R_roh_0pct, alpha=0),data.frame(R_roh_25pct, alpha=.25),
                   data.frame(R_roh_50pct, alpha=.5),data.frame(R_roh_75pct, alpha=.75),
                   data.frame(R_roh_100pct, alpha=1))

R_roh_all %>% group_by(alpha) %>% summarise_each(funs(mean))
R_roh_all %>% group_by(alpha) %>% summarise_each(funs(median))
R_roh_all %>% group_by(alpha) %>% summarise_each(funs(quantile, probs=.025))


data.frame(mean=apply(R_roh_0pct, 2, mean),
           median=apply(R_roh_0pct, 2, median),
           lowci=apply(R_roh_0pct, 2, quantile, probs=c(0.025)),
           highci=apply(R_roh_0pct, 2, quantile, probs=c(0.975)))

data.frame(mean=apply(R_roh_50pct, 2, mean),
           median=apply(R_roh_50pct, 2, median),
           lowci=apply(R_roh_50pct, 2, quantile, probs=c(0.025)),
           highci=apply(R_roh_50pct, 2, quantile, probs=c(0.975)))

data.frame(mean=apply(R_roh_100pct, 2, mean),
           median=apply(R_roh_100pct, 2, median),
           lowci=apply(R_roh_100pct, 2, quantile, probs=c(0.025)),
           highci=apply(R_roh_100pct, 2, quantile, probs=c(0.975)))





# Calculate treatment-adjusted R at specific V3s ------------------------------------------

V3 = c(0, .25, .50, .75, 1)

R_treat_0pct <- matrix(NA, nrow=nrow(fits.df), ncol=length(V3))
colnames(R_treat_0pct) <- V3
row.names(R_treat_0pct) <- 1:nrow(fits.df)
R_treat_25pct <- R_treat_50pct <- R_treat_75pct <- R_treat_100pct <- R_treat_0pct

for (i in 1:nrow(fits.df)) {
    R_treat_0pct[i,]   <- R_calc(V3=V3, Rs=fits$Rs_tmp[i], p0=fits$pk[i,1], p3=fits$pk[i,3], tau=fits$tau[i], delta=delta.2day, alpha=0)
    R_treat_25pct[i,]  <- R_calc(V3=V3, Rs=fits$Rs_tmp[i], p0=fits$pk[i,1], p3=fits$pk[i,3], tau=fits$tau[i], delta=delta.2day, alpha=.25)
    R_treat_50pct[i,]  <- R_calc(V3=V3, Rs=fits$Rs_tmp[i], p0=fits$pk[i,1], p3=fits$pk[i,3], tau=fits$tau[i], delta=delta.2day, alpha=.50)
    R_treat_75pct[i,]  <- R_calc(V3=V3, Rs=fits$Rs_tmp[i], p0=fits$pk[i,1], p3=fits$pk[i,3], tau=fits$tau[i], delta=delta.2day, alpha=.75)
    R_treat_100pct[i,] <- R_calc(V3=V3, Rs=fits$Rs_tmp[i], p0=fits$pk[i,1], p3=fits$pk[i,3], tau=fits$tau[i], delta=delta.2day, alpha=1.0)
}

R_treat_all <- rbind(data.frame(R_treat_0pct, alpha=0),data.frame(R_treat_25pct, alpha=.25),
                   data.frame(R_treat_50pct, alpha=.5),data.frame(R_treat_75pct, alpha=.75),
                   data.frame(R_treat_100pct, alpha=1))

data.frame(mean=apply(R_treat_0pct, 2, mean),
           median=apply(R_treat_0pct, 2, median),
           lowci=apply(R_treat_0pct, 2, quantile, probs=c(0.025)),
           highci=apply(R_treat_0pct, 2, quantile, probs=c(0.975)))

data.frame(mean=apply(R_treat_50pct, 2, mean),
           median=apply(R_treat_50pct, 2, median),
           lowci=apply(R_treat_50pct, 2, quantile, probs=c(0.025)),
           highci=apply(R_treat_50pct, 2, quantile, probs=c(0.975)))

data.frame(mean=apply(R_treat_100pct, 2, mean),
           median=apply(R_treat_100pct, 2, median),
           lowci=apply(R_treat_100pct, 2, quantile, probs=c(0.025)),
           highci=apply(R_treat_100pct, 2, quantile, probs=c(0.975)))


#




