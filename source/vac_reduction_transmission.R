# Vaccine reduction transmission  -----------------------------------------
if(!require('rstan')) install.packages('rstan'); library(rstan)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)

## Load in and extract data
load(file='results/R0_res_knownGT_2000iters.RData') # Loads ret_sm, fit_knGT, fit_stoch_knGT, full_data_knGT, full_data_stoch_knGT
fits <- rstan::extract(fit_knGT_fixCoV)

## Want to calculate R_{3doses}/R_{0doses} = (p3*Rs + (1-p3)*tau*Rs)/(p0*Rs + (1-p0)*tau*Rs)

# Function to calc R for 3 doses
calc.Rd <- function(p, Rs, tau){
  return(p*Rs + (1-p)*tau*Rs)
}

# Get the stan results and restrict to the 95% CI of R0
R095 <- quantile(fits$R0_preds, probs=c(0.025,0.975))
fits95 <- data.frame(R0_preds=fits$R0_preds, lRs_star=fits$lRs_star, Rs=fits$Rs_tmp, pk=fits$pk, tau=fits$tau, Vc=fits$Vc)
fits95$R0 <- fits95$pk.1*fits95$Rs + (1-fits95$pk.1)*fits95$tau*fits95$Rs
fits95 <- fits95 %>% filter(R0>=R095[1] & R0<=R095[2])

# Calculate ratio of R3d/R0d for each iteration
fits95$R3d <- NA
fits95$R0d <- NA
fits95$ratio <- NA
for (i in 1:nrow(fits95)) {
  fits95$R3d[i] <- calc.Rd(p=fits95$pk.3[i], Rs=fits95$Rs[i], tau=fits95$tau[i])
  fits95$R0d[i] <- calc.Rd(p=fits95$pk.1[i], Rs=fits95$Rs[i], tau=fits95$tau[i])
  fits95$ratio[i] <- 1- (fits95$R3d[i]/fits95$R0d[i])
  
}

mean(fits95$ratio)
quantile(fits95$ratio, probs=c(0.025, 0.975))
