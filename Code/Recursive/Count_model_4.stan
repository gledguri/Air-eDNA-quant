data {
	int time;
	int N;
	array[N] int idx;
	array[N] int<lower=0> y;
	array[N] int<lower=0> days;
}

parameters {
	array[time] real<lower=0> lambda;
}

transformed parameters{
}
model {
	for(i in 1:N){
		y[i] ~ poisson(lambda[idx[i]]*days[i]);
	}
	lambda ~ gamma(10,1);
}
