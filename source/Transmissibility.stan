data {
  int Sympt_cont; //number of contacts of symptomatic transmitters
  int Asympt_cont; //number of contacts of asymptomatic transmitters
  int Sympt_sec; //number infected by symptomatic
  int Asympt_sec; //number infected by asymptomatic
}


parameters {
  real <lower=0, upper=1> sympt_AR;
  real <lower=0, upper=1> asympt_AR;
}

transformed parameters {
  real <lower=0> tau;
  real <lower=0> X;

  X = sympt_AR/asympt_AR;
  tau = 1/(X);
  
}

model {
  
  Sympt_sec ~ binomial(Sympt_cont, sympt_AR);
  Asympt_sec ~ binomial(Asympt_cont, asympt_AR);
}
