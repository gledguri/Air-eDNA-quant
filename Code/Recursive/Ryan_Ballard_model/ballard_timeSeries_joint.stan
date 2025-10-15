data {
  int<lower=0> NqPCR; //N qPCR observations
  int<lower=0> Nt; //N time points
  int<lower=0> Nday; //N days
  int<lower=0> time_idx[NqPCR]; 
  int w[Nday]; //visual observations, standardized scale
  int<lower=0> day_idx[Nt]; //map samples onto days
  int<lower=0, upper=1> z[NqPCR]; //binary detection variable 
  vector<lower=0>[NqPCR] DNAcopies; //calculated via qPCR
  vector[Nt] meanT; //covariate: mean temperature (standardized), standardized scale
  int<lower=0> DayPhase[Nt]; //day or night
  real water_temp[Nt]; //standardized scale
  real time_to_extraction[Nt]; //standardized scale
  int<lower=0> rxn_to_day_idx[NqPCR]; //map qPCR reactions onto days
}

transformed data{
  vector[Nday] w_log; 
  
  for (d in 1:Nday){
    w_log[d] = log(w[d] + 0.5);   //transform fish counts at the locks to a log scale, adding a small value to the zeros to avoid undefined values
  }
}

parameters {
  
  // these for the true fish abundance
  real mu0; //starting point for mu in time series
  real<lower=0, upper=1> gamma; //autocorrelation slope
  vector[Nday-1] epsilon; //process variation 
  real nu; //intercept for log-scale AR1; shold be small


  // these for the molecular data
  // true DNA concentration (C) model as function of true fish (mu)
  real alpha; //intercept for true DNA concentration (log scale, should be small)
  real<lower=0> omega; //slope relating (log) fish to (log) DNA concentration on any given day
  vector[Nday] tau; //process error term relating true fish (mu, log scale) to concentration DNA (C, log scale)

  vector[4] beta; //coefficients relating detectability to true DNA concentration and to covariates
  // vector[Nt] delta; //process variability for detections
  
  real<lower=0> sigma; //observation variability on log DNA concentrations, given detection

}

transformed parameters{
  vector[Nday] mu; //true fish abundance, log scale
  vector[Nday] C; //true DNA concentration, log scale
  vector<lower=0, upper=1>[Nt] theta; //detection probability at any given time point

// true fish abundance (log scale)
  mu[1] = mu0;
  for (d in 2:Nday){
    mu[d] = nu + gamma*mu[d-1] + epsilon[d-1];
  }

// DNA
  for (d in 1:Nday){
    C[d] = alpha + omega*mu[d] + tau[d]; //mu on log scale
  }

  //observation model: detection prob; function of true concentration
  for (t in 1:Nt){
    theta[t] = inv_logit(beta[1] + beta[2]*exp(C[day_idx[t]]) + beta[3]*meanT[t] + beta[4]*DayPhase[t]); //+ delta[t]
  }

}


model {
  for (d in 1:Nday){
    w_log[d] ~ normal(mu[d], 1);  //fixed observation variability here
  }
  
  for (i in 1:NqPCR){
    z[i] ~ bernoulli(theta[time_idx[i]]);
    
    if (z[i] == 1) {
      log(DNAcopies[i]) ~ normal(C[rxn_to_day_idx[i]], sigma); 
    }
    
  }
  
  //PRIORS 
  //Fish
  mu0 ~ normal(0,5); //starting point
  nu ~ normal(0,1); //intercept
  gamma ~ beta(3,3); //autoregression
  epsilon ~ normal(0, 1); //process variation

  //DNA  
  alpha ~ normal(0, 1); //intercept
  omega ~ normal(0,5); //slope 
  tau ~ normal(0,1); //process variation

  beta[1] ~ normal(-3,1); //intercept for detectability
  beta[2] ~ gamma(1,1); //slope for detectability
  beta[3:4] ~ normal(0,5); //coeff for covariates
  // delta ~ normal(0,5); //process variation
  sigma ~ gamma(1,1); //observation variation
  
}

generated quantities{
  vector[Nday] fish_sim;
  vector[Nday] DNA_sim;
  
  for (d in 1:Nday){
    fish_sim[d] = normal_rng(mu[d], 1);
    DNA_sim[d] = normal_rng(C[d], sigma);
  }
}
