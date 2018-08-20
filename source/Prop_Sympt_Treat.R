## Want to calculate the proportion of symptomatic individuals you need to treat to get R <1. 
## R = x*(p*Rs) + (1-p)*tau*Rs 
## Where x is the proportion that remain after being treated. 
## Solve for x and find
## x = (1-((1-p)*tau*Rs))/(p*Rs)

if(!require('rstan')) install.packages('rstan'); library(rstan)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)

## Load in and extract data
load(file='results/R0_res_knownGT_2000iters.RData') # Loads ret_sm, fit_knGT, fit_stoch_knGT, full_data_knGT, full_data_stoch_knGT
fits <- rstan::extract(fit_knGT_fixCoV)


# Function to calc Vc for 3 doses
calc.x <- function(p0, Rs, tau){
    return((1-((1-p0)*Rs*tau))/(p0*Rs))
}



# For the Rohingya outbreak, we estimate containment would require --------

# Get the stan results and restrict to the 95% CI of R0
R095 <- quantile(fits$R0_preds, probs=c(0.025,0.975))
fits95 <- data.frame(R0_preds=fits$R0_preds, lRs_star=fits$lRs_star, Rs=fits$Rs[,1], pk=fits$pk, tau=fits$tau, Vc=fits$Vc)
fits95$R0 <- fits95$pk.1*fits95$Rs + (1-fits95$pk.1)*fits95$tau*fits95$Rs
fits95 <- fits95 %>% filter(R0>=R095[1] & R0<=R095[2])

# Calculate x for each iteration
fits95$x <- NA
for (i in 1:nrow(fits95)) {
    fits95$x[i] <- 1 - calc.x(p0=fits95$pk.1[i], Rs=fits95$Rs[i], tau=fits95$tau[i])
}

mean(fits95$x)
quantile(fits95$x, probs=c(0.025, 0.975))




# For all other outbreaks combined ----------------------------------------

R095 <- quantile(fits$R0_preds, probs=c(0.025,0.975))
fits95 <- data.frame(R0_preds=fits$R0_preds, lRs_star=fits$lRs_star, Rs=fits$Rs_tmp, pk=fits$pk, tau=fits$tau, Vc=fits$Vc)
fits95$R0 <- fits95$pk.1*fits95$Rs + (1-fits95$pk.1)*fits95$tau*fits95$Rs
fits95 <- fits95 %>% filter(R0>=R095[1] & R0<=R095[2])

# Calculate x for each iteration
fits95$x <- NA
for (i in 1:nrow(fits95)) {
    fits95$x[i] <- 1- calc.x(p0=fits95$pk.1[i], Rs=fits95$Rs[i], tau=fits95$tau[i])
}

mean(fits95$x)
quantile(fits95$x, probs=c(0.025, 0.975))
