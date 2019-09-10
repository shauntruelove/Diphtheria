
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('rstan')) install.packages('rstan'); library(rstan)
if(!require('lubridate')) install.packages('lubridate'); library(lubridate)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('scales')) install.packages('scales'); library(scales)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # Detect and set to multi cores
source('source/restrict_stan_input_data.R')


# DATA --------------------------------------------------------------------

# Re-run data source
data.rerun <- FALSE
if (data.rerun) {
  source('./source/R0_outbreak_data.R')
  STAN_data <- get_STAN_data()
}

load(file='./data/R0_STAN_data.RData') #Loads full_data
load(file='./data/R0_STAN_data_stoch.RData') #Loads full_data_stoch
load(file="./results/gentime_params.RData")  #Loads gentime_params

# Add additional variables to the STAN list object
length(full_data)
full_data[['shape_mean']] <- gentime_params$shape_mean
full_data[['shape_sd']] <- gentime_params$shape_sd
full_data[['scale_mean']] <- gentime_params$scale_mean
full_data[['scale_sd']] <- gentime_params$scale_sd
full_data[['logtau_mean']] <- -1.437 
full_data[['logtau_sd']] <- 0.09 #0.285
full_data[['p0_mean']] <- .69
full_data[['p0_sd']] <- 0.097
full_data[['VE_mean']] <- c(.705, .872) 
full_data[['VE_sd']] <- c(.110, .047) 
full_data[["Vk_cov"]] <- c(0.10,0.10,0.10)

full_data_stoch[['shape_mean']] <- gentime_params$shape_mean
full_data_stoch[['shape_sd']] <- gentime_params$shape_sd
full_data_stoch[['scale_mean']] <- gentime_params$scale_mean
full_data_stoch[['scale_sd']] <- gentime_params$scale_sd
full_data_stoch[['logtau_mean']] <- -1.437 
full_data_stoch[['logtau_sd']] <- 0.09 #0.285
full_data_stoch[['p0_mean']] <- .69
full_data_stoch[['p0_sd']] <- 0.097
full_data_stoch[['VE_mean']] <- c(.705, .872) 
full_data_stoch[['VE_sd']] <- c(.110, .047) 
full_data_stoch[["Vk_cov"]] <- c(0.10,0.10,0.10)

# Fix NAs in N_mat
full_data$N_mat[is.na(full_data$N_mat)] <- 0
full_data_stoch$N_mat[is.na(full_data_stoch$N_mat)] <- 0
full_data_knGT <- full_data
full_data_stoch_knGT <- full_data_stoch

full_data_knGT_noRoh <- restrict_stan_input_data(data=full_data_knGT, outbreaks_remove=c(1))


# MODEL ---------------------------------------------------------------

# Known GT model - Fixed CoV
fileName_knownGT_fixedCoV <- "./source/R0stan_knownGT_fixedCoV.stan"
ret_known <- stanc(fileName_knownGT_fixedCoV) # Check Stan file
ret_sm_knownGT <- stan_model(stanc_ret = ret_known) # Compile Stan code
save(fileName_knownGT_fixedCoV, ret_known, ret_sm_knownGT, file="./source/R0stan_knownGT_fixedCoV_compiled.RData")
load(file="./source/R0stan_knownGT_fixedCoV_compiled.RData")


# RUN TEST ---------------------------------------------------------------

fit_knGT <- sampling(ret_sm_knownGT, warmup=100, iter=200, seed=123, data=full_data_knGT, chains=4, control=list(adapt_delta=.85))
fit_knGT

fit_knGT_noRoh <- sampling(ret_sm_knownGT, warmup=100, iter=200, seed=123, data=full_data_knGT_noRoh, chains=4, control=list(adapt_delta=.85))
fit_knGT_noRoh

# fit_stoch_knGT <- sampling(ret_sm_knownGT, warmup=100, iter=200, seed=123, data=full_data_stoch_knGT, chains=4, control=list(adapt_delta=.85))
# fit_stoch_knGT

# fit <- sampling(ret_sm, warmup=100, iter=200, seed=123, data=full_data, chains=4, control=list(adapt_delta=.85))
# fit
# 
# fit_stoch <- sampling(ret_sm, warmup=100, iter=200, seed=123, data=full_data_stoch, chains=4, control=list(adapt_delta=.85))
# fit_stoch


# Save STAN results
#save(ret_sm, fit, fit_stoch, full_data, full_data_stoch, file='./trunk/Results/R0res_200iters.RData')
#save(ret_sm_knownGT, fit_knGT, fit_stoch_knGT, full_data_knGT, full_data_stoch_knGT, file='./trunk/Results/R0res_knownGT_200iters.RData')
save(ret_sm_knownGT, fit_knGT, full_data_knGT, fit_knGT_noRoh, full_data_knGT_noRoh, file='./results/R0res_knownGT_fixedCoV10prct_200iters.RData')
load(file='./results/R0res_knownGT_fixedCoV10prct_200iters.RData') # Loads ret_sm, fit, full_data, full_data_stoch
fit_knGT_fixCoV <- fit_knGT
load(file='./results/R0res_knownGT_200iters.RData') # Loads ret_sm, fit, full_data, full_data_stoch


fit_knGT_fixCoV
fit_knGT


#R0.posterior <- extract(fit)$R0  
summary(exp(extract(fit_knGT)$lRs_star))
mean(extract(fit_knGT)$R)


fit_knGT


# RUN FULL ---------------------------------------------------------------
warmup <- 1000
iter <- 2000

fit_knGT <- sampling(ret_sm_knownGT, warmup=warmup, iter=iter, seed=123, data=full_data_knGT, chains=4, control=list(adapt_delta=.95))
fit_knGT

# fit_stoch_knGT <- sampling(ret_sm_knownGT, warmup=warmup, iter=iter, seed=123, data=full_data_stoch_knGT, chains=4, control=list(adapt_delta=.95))
# fit_stoch_knGT
# 
# # Unknown GT
# fit <- sampling(ret_sm, warmup=warmup, iter=iter, seed=123, data=full_data, chains=4, control=list(adapt_delta=.95))
# fit
# 
# fit_stoch <- sampling(ret_sm, warmup=warmup, iter=iter, seed=123, data=full_data_stoch, chains=4, control=list(adapt_delta=.95))
# fit_stoch


# Save STAN results
# save(ret_sm, fit, fit_stoch, full_data, full_data_stoch, file='./results/R0_res_1000iters.RData')
# save(ret_sm, fit_knGT, fit_stoch_knGT, full_data_knGT, full_data_stoch_knGT, file='./results/R0res_knownGT_2000iters.RData')
save(ret_sm, fit_knGT, full_data_knGT, file='./results/R0res_knownGT_fixedCoV10prct_2000iters.RData')
load(file='./results/R0res_knownGT_fixedCoV10prct_2000iters.RData') # Loads ret_sm, fit_knGT, fit_stoch_knGT, full_data, full_data_stoch





# ANALYSIS ----------------------------------------------------------------

chains <- extract(fit_knGT)

plot(chains$R0_preds, chains$Vc, ylim=c(-1,2))
abline(h=1, col='red',lty=3)
abline(h=0, col='grey',lty=1)











fit <- fit_knGT

chains <- extract(fit_knGT)

summary(chains$Vc)
summary(chains$tau)



R0ests <- rstan::extract(fit_knGT)$R0_preds
Rfits <- rstan::extract(fit_knGT)$R



# Compare vaccination coverages priors and posterior distributions

fits <- rstan::extract(fit_knGT)

V_est <- fits$V_est
V_input <- full_data_knGT$Vk_pop

for (i in 1:22){
  V_tmp <- V_est[,i,]
  plot(density(V_tmp[,1]), col='red', lwd=2, main=paste0("Outbreak-",i), xlim=c(0,1), ylim=c(0,30))
  lines(density(V_tmp[,2]), col='blue', lwd=2)
  lines(density(V_tmp[,3]), col='darkgreen',lwd=2)
  abline(v=V_input[i,1], col='red', lty=3)
  abline(v=V_input[i,2], col='blue', lty=3)
  abline(v=V_input[i,3], col='darkgreen', lty=3)
  text(.2,25,paste0("V0=",V_input[i,1]))
}


plot(fits$R0_preds, fits$Vc)
abline(h=1, lty=3)

mean(V_est[,,1])









pk_ests <- rstan::extract(fit_knGT)$pk
p3ests <- pk_ests[,3]
summary(p3ests)
hist(p3ests)
1-summary(p3ests)
1-quantile(p3ests, probs=c(0.025, 0.975))

1-((1-.872)*0.71)

Vc <- rstan::extract(fit_knGT)$Vc
summary(Vc)
quantile(Vc, probs=c(0.025, 0.975))





Rs_tmp = R0_preds / (p0 + (1-p0)*tau)
Vc = (p0 + tau - tau*p0 - 1/Rs_tmp) / (p0 - tau*p0 - p3 + p3*tau); 
Vc = (p0 + tau - tau*p0 - (p0 + (1-p0)*tau)/R0_preds) / (p0 - tau*p0 - p3 + p3*tau)
Vc = (p0 + tau - tau*p0 - ((p0*(R0_preds / (p0*Rs_tmp + (1-p0)*tau*Rs_tmp) + (1-p0)*tau*Rs_tmp) / R0_preds) / (p0 - tau*p0 - p3 + p3*tau)))








plot(fit)
fit
fits_ <- extract(fit)

# R0 Estimates
summary(exp(extract(fit)$lRstar))
exp(quantile(extract(fit)$lRstar, probs = c(.025, .975)))

mean(fits_$R)

# R0 Predictions
R0.posterior <- extract(fit)$R0_preds  
summary(R0.posterior)
quantile(R0.posterior, probs = c(.025, .975))


# R effective
R_effs <- extract(fit)$R
summary(R_effs)
summary(as.numeric(R_effs))
range(as.numeric(R_effs))
mean_Rs <- colMeans(R_effs)
range(mean_Rs)
quantile(as.numeric(R_effs), probs = c(.025, .975))
mean(mean_Rs)
quantile(as.numeric(mean_Rs), probs = c(.025, .975))


calc_R0 <- function(Rs, p0=0.8) {
  p0*Rs + (1-p0)*.28*Rs
}



# Extract R0 estimates from each outbreak
R0_ests <- extract(fit)$R / extract(fit)$R0_adj 

hist(extract(fit)$R[,1])
hist(R0_ests[,1])
hist(R0_ests[,2])
hist(R0_ests[,3])

summary(R0_ests)



# Check VE

hist(extract(fit)$VE)



