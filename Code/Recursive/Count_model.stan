data {
  int time;
  array[time] int<lower=0> y;
}

parameters {
	array[time] real<lower=0> lambda;
	// real<lower=0> lambda;
}

transformed parameters{
	array[time] real lambda_s;
	lambda_s[1] = (lambda[1]+lambda[2])/2;
	lambda_s[time] = (lambda[time]+lambda[(time-1)])/2;
		for (t in 2:(time-1)){
			lambda_s[t] = (lambda[t-1] + lambda[t] + lambda[t+1])/3;
		}
}
model {
		// target += poisson_cdf(y|lambda);
		y ~ poisson(lambda_s);
		lambda ~ gamma(10,1);
}

generated quantities {
  array[time] int y_pred;
  array[time] int y_step;
  for (t in 1:time) {
    y_step[t] = poisson_rng(lambda_s[t]);  // Generate counts at each time step (single value)
    y_pred[t] = (t == 1) ? y_step[t] : y_pred[t - 1] + y_step[t]; // Cumulative sum
  }
}
