
data {
  int<lower=1> N;               // Number of subjects
  int<lower=1> J;               // the number of groups (Number of studies from which data comes)
  int<lower=1,upper=J> id[N];   // vector of group ids (indices)
  int<lower=0> K;               // Number of predictor parameters (fixed)
  matrix[N,K] X;                // fixed effects model integer array
  int<lower=0,upper=1> Y[N];    // response variable - Outcome for each subect (died=1, survived=0)
}

transformed data {
  int<lower=0> K_real;
  K_real = K - 1;
}

parameters {
  real beta0;
  vector[K_real] beta;               // Fixed effects
  real<lower=0, upper=100> sigma;    // SDs for random effect
  vector[J] beta_j;                  // Study Intercept (CFR coefficient in each study)
}

model {
  // Priors on group parameters
  beta0 ~ normal(0, 100);
  beta ~ normal(0, 100);

  beta_j ~ normal(beta0, sigma);
  
  // The likelihood 
  Y ~ bernoulli_logit(beta_j[id] + beta[1]*X[,2]);

}

generated quantities {         // simulate quantities of interest
  real CFR0;
  real CFR1;

  CFR0 = inv_logit(beta0);           
  CFR1 = inv_logit(beta0 + beta[1]); 
}

