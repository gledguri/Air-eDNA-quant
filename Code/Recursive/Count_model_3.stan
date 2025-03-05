data {
	int time;
	// int N;
	// array[N] int idx;
	array[time] int<lower=0> y;
	array[time] int<lower=0> days;
	array[2] real lambda_rp;
	array[2] real lambda_p;
	array[2] real beta_p;
	array[2] real epsilon_p;
}

parameters {
	array[time] real<lower=0> lambda_raw;
	real<lower=0> lambda_0;
	real<lower=0, upper=1> beta;
	array[time-1] real epsilon;
}

transformed parameters{
  array[time] real<lower=0> lambda;
  lambda[1] = lambda_0;

  for (t in 2:time){
    // lambda[t] = (beta*lambda[t-1] + epsilon[t-1]);
    lambda[t] = exp(log(beta*lambda[t-1]) + epsilon[t-1]);
  }
}

model {
	for(i in 1:time){
		y[i] ~ poisson(lambda[i]*days[i]);
		y[i] ~ poisson(lambda_raw[i]*days[i]);
	}
	// Priors
	// lambda ~ gamma(10,1);
	lambda_raw ~ gamma(lambda_rp[1],lambda_rp[2]);
	lambda_0 ~ gamma(lambda_p[1],lambda_p[2]);
  beta ~ beta(beta_p[1],beta_p[2]); //autoregression
  epsilon ~ normal(epsilon_p[1], epsilon_p[2]); //process variation
}
