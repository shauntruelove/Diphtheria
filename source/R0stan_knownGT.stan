
data {
  // Data inputs
  int M;                                    // number of outbreaks
  int<lower=1> T[M];                        // maximum time of each outbreak/exponential growth curve
  int<lower=1> T_max;                       // maximum length among all outbreaks
  int N_mat[M,T_max];                       // all cases at each time step, for each outbreak, minus t=0
  int<lower=1> N0[M];                       // Initial cases
  matrix[M,3] Vk_pop;                       // proportion receiving k doses in each outbreak
  int<lower=0, upper=1> use_schick[M];      // Logical, if outbreak occurred during prevaccination era (pre-1923)
  int<lower=0, upper=1> est_Vs[M];          // If stan should be allowed to estimate Vs for population
  int<lower=0, upper=1> alpha[M];           // If cases were treated with antibiotics in the outbreak (logical 0/1)

  // Parameter inputs: Generation Time
  real<lower=0> shape_mean;           // Shape parameter for Generation Time gamma distribution [MEAN]
  real<lower=0> scale_mean;           // Scale parameter for Generation Time gamma distribution [MEAN]

  // Parameter inputs: R0 and Transmission
  real<upper=-0.001> logtau_mean;           // Ratio of R0 in asymptomatics vs syptomatics
  real<lower=0> logtau_sd;                  // Ratio of R0 in asymptomatics vs syptomatics
  real<lower=0, upper=1> p0_mean;           // Mean prop symptomatic among unvaccinated.
  real<lower=0, upper=1> p0_sd;             // sd prop symptomatic among unvaccinated.
  real<lower=0, upper=1> VE_mean[2];        // Mean VE with k doses.
  real<lower=0, upper=1> VE_sd[2];          // SD VE with k doses.
  real<lower=0, upper=1> delta_range[2];    // range of reduction ratio due to treatment.

}

transformed data {
  int<lower=0> T_m1[M];         // Max time T for each outbreak, minus 1. This is done as j<T in White Pagano model
  int<lower=0> T_max_m1;   	    // Max time T among all outbreaks, minus 1. This is done as j<T in White Pagano model
  real<lower=0> Vk_cov[3];      // Coefficient of variation of the inputted Vs
  real<lower=0.0001> V3_sd[M];   // Normalized SD for V3 for each outbreak (have to explicitly generate it because it cannot = 0)
  real<lower=0.0001> V1_sd[M];   // Normalized SD for V1 for each outbreak (have to explicitly generate it because it cannot = 0)

  for (m in 1:M) { 
    T_m1[m] = T[m]-1;                       // Generate T_m1 for each outbreak
  }
  T_max_m1 = T_max-1;                       // Generate T_max_m1
  
  // Coefficient of Variation - Vk
  for (i in 1:3){
    Vk_cov[i] = sd(Vk_pop[,i]) / mean(Vk_pop[,i]);
  }
  // Because SD cannot be 0 for a normal distrib, we are adding 0.0001 to make it work
  for (m in 1:M){
    V3_sd[m] = Vk_pop[m,3]*Vk_cov[3] + 0.0001;
    V1_sd[m] = Vk_pop[m,1]*Vk_cov[1] + 0.0001;
    
    // If the true vaccination coverage is known
    if (est_Vs[m]==0){
      V1_sd[m] = 0.0001;
      V3_sd[m] = 0.0001;
    }
  }
    
}


parameters {
  // R0 Parameters
  real lRs[M];                              // log Rs from each outbreak
  real<lower=0> lRs_star;                   // log Rs (mean accross all oubreaks) [FIT]
  real<lower=0> sigma;                      // Standard deviation in log R0 accross outbreaks [FIT]

  // Diphtheria dynamics Parameters
  real<upper=-0.001> logtau;                // Ratio of R0 in asymptomatics vs syptomatics
  real<lower=0.001,upper=1> VE_1;
  real<lower=0.001,upper=1> VE_3;
  real<lower=0.001,upper=1> p0;             // Sampled probability of symptoms - unvaccinated
  real<lower=0,upper=1> delta;              // relative transmissibility per treated symptomatic vs untreated
  simplex[3] V_est[M];                      // estimated proportion receiving 0 doses in each outbreak
}   


transformed parameters {
  real<lower=0> R[M];                       // Effective R from outbreaks
  real<lower=0> Rs[M];                      // Reproductive number among symptomatic cases
  real<lower=0,upper=1> pk[3];              // Prop symptomatic after receiving k doses of vaccine.
  real<lower=0.001,upper=1> tau;            // Ratio of R0 in asymptomatics vs syptomatics

  tau = exp(logtau);

  // calculate the proportion symptomatic from VE estimates for 1-2 and 3+ doses
  pk[1] = p0;
  pk[2] = (1-VE_1) * pk[1];
  pk[3] = (1-VE_3) * pk[1];

  // exponentiate lRs
  for (m in 1:M) {
    Rs[m] = exp(lRs[m]);                   
  }

  // Calculate Rs from R for each outbreak from sampled parameters
  for (m in 1:M) { 
    R[m] = 0;
    if (use_schick[m]==1){
      R[m] = Rs[m]*(V_est[m,3]*(1 - alpha[m] + alpha[m]*delta) + (1-V_est[m,3])*tau);
    } else {
      R[m] = Rs[m]*(V_est[m,1]*(pk[1]*(1 - alpha[m] + alpha[m]*delta) + (1-pk[1])*tau) + 
                    V_est[m,2]*(pk[2]*(1 - alpha[m] + alpha[m]*delta) + (1-pk[2])*tau) + 
                    V_est[m,3]*(pk[3]*(1 - alpha[m] + alpha[m]*delta) + (1-pk[3])*tau));
    }
  }
  
}


model { 
  matrix[M,T_max_m1] p;                        // Generation time probability distribution
  real p_sum;
  real mu_tmp;

  // Priors (weakly informative)
  lRs_star ~ normal(0, 100);                
  sigma ~ normal(0,100);                       
  p = rep_matrix(0, M, T_max_m1);               // generation time prob distribution
  delta ~ uniform(delta_range[1], delta_range[2]);

  // Estimate global log Rs_star, with error (sigma), from log Rs
  lRs ~ normal(lRs_star, sigma);
  
  // Estimate the Vaccination Coverages for each population
  for (m in 1:M){
    V_est[m,3] ~ normal(Vk_pop[m,3], V3_sd[m]);
    V_est[m,1] ~ normal(Vk_pop[m,1], V1_sd[m]);
  }

  // R0 Adjustment parameters
  logtau ~ normal(logtau_mean, logtau_sd);
  target += normal_lpdf(VE_1 | VE_mean[1], VE_sd[1]);
  target += normal_lpdf(VE_3 | VE_mean[2], VE_sd[2]);
  target += normal_lpdf(p0 | p0_mean, p0_sd);


   // Calculate p vector -- the generation time probability distribution
     for (m in 1:M) {
        for (j in 1:T_m1[m]) {
          p[m,j] = gamma_cdf(j, shape_mean, 1/scale_mean) - gamma_cdf(j-1, shape_mean, 1/scale_mean);  
        }
     }

    // Normalize p vector to sum to 1 (because it's truncated in some cases - per White and Pagano pp.3005)
     for (m in 1:M) {
        p_sum = sum(p[m,1:T_max_m1]);
        for (j in 1:T_max_m1) {
          p[m,j] = p[m,j] / p_sum;
        }
     }
	
  // Calculate mu for each outbreak
  for (m in 1:M) {
      
    // Create mu vector here because each outbreak has different T
    vector[T[m]] mu;                             
    
    for (t in 1:T[m]) {
      mu_tmp = 0.0;
      for (j in 1:min(T_m1[m], t)) {
  	    if (t!=j) {
  	        mu_tmp = N_mat[m][t-j]*p[m,j] + mu_tmp;  
  	    } else {
  	        mu_tmp = N0[m]*p[m,j] + mu_tmp; // If t==j, we use N0 
  	    }
      }
      mu[t] = mu_tmp*R[m]; // Multiply mu vector by outbreak-specific R
    }
    	
    // likelihood of incidence vector, given vector of mu
    N_mat[m,1:T[m]] ~ poisson(mu);             
  }
}


generated quantities {
    // create random R0
    real R0_preds;
    real Rs_tmp;
    real Vc;
    
    Rs_tmp = exp(normal_rng(lRs_star, sigma));
    
    R0_preds = pk[1]*Rs_tmp + (1-pk[1])*tau*Rs_tmp;
    
    Vc = (pk[1] + tau - tau*pk[1] - 1/Rs_tmp) / (pk[1] - tau*pk[1] - pk[3] + pk[3]*tau); 
}
