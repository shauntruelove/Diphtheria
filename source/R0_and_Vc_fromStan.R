
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('rstan')) install.packages('rstan'); library(rstan)
if(!require('lubridate')) install.packages('lubridate'); library(lubridate)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('scales')) install.packages('scales'); library(scales)

source('source/restrict_stan_input_data.R')
options(max.print = .Machine$integer.max)


# LOAD STAN RESULTS --------------------------------------------------------------------

load(file='results/R0_res_knownGT_2000iters.RData') # Loads ret_sm_knownGT, fit_knGT, fit_knGT_fixCoV, fit_knGT_fixCoV_noRoh, full_data, full_data_knGT


# View Results ------------------------------------------------------------

fit_knGT_fixCoV  # Final 
fit_knGT



# Get Results for Paper ---------------------------------------------------

fits <- rstan::extract(fit_knGT_fixCoV)

# R0 Estimates
summary(fits$R0_preds)
quantile(fits$R0_preds, probs=c(0.025,.975))


# Vc Estimates
summary(fits$Vc)
quantile(fits$Vc, probs=c(0.025,.975))


# R Estimates - 23 Outbreaks
R <- fits$R
colMeans(R)
apply(R, 2, median)
range(apply(R, 2, median))
mean(apply(R, 2, median))
summary(apply(R, 2, median))
quantile(apply(R, 2, median), probs=c(0.025,.975))










# Compare Results to the results with weak priors for vaccination coverage ------------------------------------------------------------

fit_knGT_fixCoV  # Final keeper results
fit_knGT


















