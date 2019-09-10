
data {
  int<lower=1> N;               // Number of subjects
  int<lower=1> J;               // the number of groups (Number of studies from which data comes)
  int<lower=1,upper=J> id[N];   // vector of group ids (indices)
  int<lower=0> K;               // Number of predictor parameters (fixed)
  matrix[N,K] X;                // fixed effects model integer array
  int<lower=0,upper=1> Y[N];    // response variable - Outcome for each subect (died=1, survived=0)
}

parameters {
  real beta;                         // Fixed effect
  real<lower=0, upper=100> sigma;    // SDs for random effect
  vector[J] beta_j;                  // Study Intercept (VE coefficient in each study)
  vector[J] beta0;                   // Study Intercept (VE coefficient in each study)

}

model {
  // Priors on group parameters
  beta0 ~ normal(0, 10);
  beta ~ normal(0, 10);
  
  beta_j ~ normal(beta, sigma);

  // The likelihood 
  for (n in 1:N) {
    Y[n] ~ bernoulli_logit(beta0[id[n]] + beta_j[id[n]]*X[n,2]);
  }
}

generated quantities {         // simulate VE
  real VE;
  VE = 1 - exp(beta); 
}
