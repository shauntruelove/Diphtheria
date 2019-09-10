
data {
  int<lower=1> N;               // Number of subjects
  int<lower=1> J;               // the number of groups (Number of studies from which data comes)
  int<lower=1,upper=J> id[N];   // vector of group ids (indices)
  int<lower=0> K;               // Number of predictor parameters (fixed)
  //int<lower=0,upper=1> X[N,K];  // fixed effects model integer array
  matrix[N,K] X;                // fixed effects model integer array
  int<lower=0,upper=1> Y[N];    // response variable - Outcome for each subect (died=1, survived=0)
}

transformed data {
  int<lower=0> K_real;
  // matrix[N,K] X_mat;
  // matrix[N,J] Z_mat;
  K_real = K - 1;
  // X_mat = to_matrix(X);
  // Z_mat = to_matrix(Z);
}

parameters {
  real beta0;
  vector[K_real] beta;               // Fixed effects
  real<lower=0, upper=100> sigma;    // SDs for random effect
  vector[J] beta_j;                  // Study Intercept (CFR coefficient in each study)
}

// transformed parameters {
// 
//   for (i in 1:N) {
//     mu[i] = beta_j[id[[i]]] + beta[1]*X[i,2] + beta[2]*X[i,3];
//   }
// }

model {
  //real mu[N];    

  // Priors on group parameters
  beta0 ~ normal(0, 100);
  beta ~ normal(0, 100);

  beta_j ~ normal(beta0, sigma);
  
  // The likelihood 
  // for (i in 1:N) {
  //   Y ~ bernoulli_logit(beta_j[id[i]] + beta[1]*X[i,2] + beta[2]*X[i,3]);
  // }
  Y ~ bernoulli_logit(beta_j[id] + beta[1]*X[,2] + beta[2]*X[,3]);

}

generated quantities {         // simulate quantities of interest
  real CFR0;
  real CFR1;
  real CFR2;

  CFR0 = inv_logit(beta0);           
  CFR1 = inv_logit(beta0 + beta[1]); 
  CFR2 = inv_logit(beta0 + beta[2]); 

}

