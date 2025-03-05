data {
	int time;
	int N;
	array[N] int idx;
	array[N] int<lower=0> y;
}

parameters {
	array[time] real<lower=0> lambda;
	// real<lower=0> lambda;
}

transformed parameters{
}
model {
	y ~ poisson(lambda[idx]);
	lambda ~ gamma(10,1);
}
