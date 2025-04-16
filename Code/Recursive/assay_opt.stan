data {
  int N_st_q; // Total number of observations in qPCR standard samples
  int N_st_qp; // Total number of observations in qPCR standard samples for only detected samples
  array[N_st_q] int Z_qst; // Presence/Absence response of qPCR standard data
  array[N_st_q] real S_q; // Known concentration (log10) in qPCR data
  array[N_st_qp] real R_qst; // Ct values of qPCR standard data for only detected samples
  array[N_st_qp] real S_q_p; // Known concentration (log10) in qPCR data for only detected samples
}
parameters {
  real alpha_0;
  real alpha_1;
  real eta_0;
  real eta_1;
  real gamma_0;
  real<upper=0> gamma_1;
}
transformed parameters {
  vector[N_st_q] theta_st;
  vector[N_st_qp] mu_st;
  vector[N_st_qp] sigma_st;

  for (i in 1:N_st_q) {
    theta_st[i] = alpha_0 + (alpha_1 * (S_q[i] - 3));
  }
  for (i in 1:N_st_qp) {
    mu_st[i] = eta_0 + (eta_1 * S_q_p[i]);
    sigma_st[i] = exp(gamma_0 + (gamma_1 * S_q_p[i]));
  }
}
model {
	// Weighted contributions to the posterior
	Z_qst ~ bernoulli(inv_logit(theta_st));
  R_qst ~ normal(mu_st,sigma_st);
  // target += (1-inv_logit(theta_st)) * bernoulli_lpmf(Z_qst | inv_logit(theta_st));
  // target += 1 * normal_lpdf(R_qst | mu_st, sigma_st);
  print("Log-posterior: ", target());
  // Priors
  alpha_0 ~ normal(0, 1);
  alpha_1 ~ normal(0, 1);
  eta_0 ~ normal(0, 10);
  eta_1 ~ normal(0, 3);
  gamma_0 ~ normal(1, 0.1);
  gamma_1 ~ normal(0, 0.1);
}
