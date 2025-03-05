functions {
  array[] real average_by_idx(int N, int K, array[] real input, array[] int index) {
    array[K] real output;
    array[K] int count;
    count = rep_array(0, K);
    output = rep_array(0.0, K);
    for (i in 1:N) {
      count[index[i]] = count[index[i]] + 1;
      output[index[i]] = output[index[i]] + input[i];
    }
    for (i in 1:K) {
      if (count[i] != 0) { // Avoid division by zero
        output[i] = output[i] / count[i];
      }
    }
    return output;
  }
}
// 
data {
// 	//Numbers of dimentions
	int N_st_q; // Total number of observation in qPCR standard samples
	int N_en_wat_q; // Total number of observation in qPCR environmental samples
	int N_en_air_q; // Total number of observation in qPCR environmental samples
// 	// 
	int N_st_qp; // Total number of observation in qPCR standard samples for only detected samples
	int N_en_wat_qp; // Total number of observation in qPCR environmental samples for only detected samples
	int N_en_air_qp; // Total number of observation in qPCR environmental samples for only detected samples
// 	// 
  int N_filt; // Total number of plates used
  int N_plate; // Total number of plates used
	int N_j_wat; // Number of samples in environmental data
	int N_j_air; // Number of samples in environmental data
	int N_ft_air; //
// 	//
// 	// // Indexes
//   // // // Plate index
  array[N_st_qp] int plate_st_idx; // Plate index
  array[N_en_wat_qp] int plate_en_wat_idx; // Plate index for water samples
  array[N_en_air_qp] int plate_en_air_idx; // Plate index for air samples
// 	// // // Binomial model
	array[N_en_wat_q] int j_qen_wat_idx; // Species and standard index for qPCR environmental samples
	array[N_en_air_q] int j_qen_air_idx; // Species and standard index for qPCR environmental samples
// 	// // // Continious model
	array[N_en_wat_qp] int j_qen_wat_p_idx; // Species and standard index for qPCR environmental samples
	array[N_en_air_qp] int j_qen_air_p_idx; // Species and standard index for qPCR environmental samples
// 	// Data
// 	// // // Binomial model
	array[N_st_q] int Z_qst; // Presence/Absence response of qPCR standard data
	array[N_en_wat_q] int Z_qen_wat; // Presence/Absence response of qPCR environmental data
	array[N_en_air_q] int Z_qen_air; // Presence/Absence response of qPCR environmental data
// 	// // // Continious model
	array[N_st_qp] real R_qst; // Ct values of qPCR standard data for only detected samples
	array[N_en_wat_qp] real R_qen_wat; // Ct values of qPCR environmental data for only detected samples
	array[N_en_air_qp] real R_qen_air; // Ct values of qPCR environmental data for only detected samples
// 	// // // 
	array[N_st_q] real S_q; // Known concentration (log10) in qPCR data
	array[N_st_qp] real S_q_p; // Known concentration (log10) in qPCR data for only detected samples
// 	// 
	int time;
	array[time] int<lower=0> N;
	array[time] int<lower=0> E;
	array[N_j_wat] int t1_idx;
// 	#array[N_j_wat] int t0_idx;
	// array[N_j_air] int air_wat_idx;
// 	array[N_j_air] int air_filt_idx;
	array[N_j_air] int ft_idx;
	array[N_ft_air] int ft_idx2;
	array[N_ft_air] int air_wat_idx2;
	array[N_ft_air] int air_filt_idx2;
	real tau_p1;
	real tau_p2;
	real logit_phi_mu;
	real logit_phi_sd;
}
// 
parameters {
	////////////////////////////////////////////////////////////////////// Bernuli parameters
	real logit_phi;
	////////////////////////////////////////////////////////////////////// Continous parameters
  vector[N_plate] beta_0;
	real beta_1;
	real gamma_0;
	real<upper=0> gamma_1;
	////////////////////////////////////////////////////////////////////// Count parameters
	// array[time] real<lower=0> lambda;
	////////////////////////////////////////////////////////////////////// Count to Water param
	real omega;
	// ////////////////////////////////////////////////////////////////////// Water to Air param
	array[N_filt] real tau;
	array[N_ft_air] real epsilon;
	array[N_filt] real eta;
	vector[N_j_air] log_A;
	// vector[N_j_wat] log_W;
}
// 
transformed parameters{
	////////////////////////////////////////////////////////////////////// Bernuli Model
	//////////////////////////////////////////////////////////////////////// Standards
	real phi;
	phi = inv_logit(logit_phi);
  vector[N_st_q] p_tmp_st;
  vector[N_st_q] psi_st;
  // 
	for (i in 1:N_st_q){
		p_tmp_st[i] = exp(S_q[i]) * -phi;
    psi_st[i] = log1m_exp(p_tmp_st[i]) - p_tmp_st[i];
	}
	//////////////////////////////////////////////////////////////////////// Air samples
  vector[N_en_air_q] p_tmp_un_air;
  vector[N_en_air_q] psi_un_air;
  // 
	for (i in 1:N_en_air_q){
		p_tmp_un_air[i] =  exp(log_A[j_qen_air_idx[i]]) * -phi;
    psi_un_air[i] = log1m_exp(p_tmp_un_air[i]) - p_tmp_un_air[i];
	}
	////////////////////////////////////////////////////////////////////// Continuous Model
	//////////////////////////////////////////////////////////////////////// Standards
	vector[N_st_qp] mu_st;
	vector[N_st_qp] sigma_st;
	// 
	for (i in 1:N_st_qp){
		mu_st[i] = beta_0[plate_st_idx[i]] + (beta_1 * S_q_p[i]);
		sigma_st[i] = exp(gamma_0+(gamma_1 * S_q_p[i]));
	}
	// //////////////////////////////////////////////////////////////////////// Air samples
	vector[N_en_air_qp] mu_en_air;
	vector[N_en_air_qp] sigma_en_air;
	// 
	for (i in 1:N_en_air_qp){
		mu_en_air[i] = beta_0[plate_en_air_idx[i]] + (beta_1 * log_A[j_qen_air_p_idx[i]]);
		sigma_en_air[i] = exp(gamma_0+(gamma_1 * log_A[j_qen_air_p_idx[i]]));
	}
	//////////////////////////////////////////////////////////////////////// Mean Air
	vector[N_ft_air] mean_log_A;
	vector[N_j_air] delta;
	// 
	mean_log_A = log(to_vector(average_by_idx(N_j_air, N_ft_air, to_array_1d((exp(log_A))), ft_idx)));
	for (i in 1:N_j_air){
		delta[i] = mean_log_A[ft_idx[i]] - log_A[i];
	}
	//////////////////////////////////////////////////////////////////////// Air to Water
	vector[N_j_wat] log_W;
	// 
	for (i in 1:N_ft_air){
		log_W[air_wat_idx2[i]] = eta[air_filt_idx2[i]] + mean_log_A[ft_idx2[i]] + epsilon[ft_idx2[i]];
	}
	//////////////////////////////////////////////////////////////////////// Water samples
  vector[N_en_wat_q] p_tmp_un_wat;
  vector[N_en_wat_q] psi_un_wat;
  // 
	for (i in 1:N_en_wat_q){
		p_tmp_un_wat[i] =  exp(log_W[j_qen_wat_idx[i]]) * -phi;
    psi_un_wat[i] = log1m_exp(p_tmp_un_wat[i]) - p_tmp_un_wat[i];
	}
	// //////////////////////////////////////////////////////////////////////// Water samples
	vector[N_en_wat_qp] mu_en_wat;
	vector[N_en_wat_qp] sigma_en_wat;
	// 
	for (i in 1:N_en_wat_qp){
		mu_en_wat[i] = beta_0[plate_en_wat_idx[i]] + (beta_1 * log_W[j_qen_wat_p_idx[i]]);
		sigma_en_wat[i] = exp(gamma_0+(gamma_1 * log_W[j_qen_wat_p_idx[i]]));
	}
	// ////////////////////////////////////////////////////////////////////// Count to Water
	array[time] real <lower=0> lambda;
	// 
	for (i in 1:N_j_wat){
		lambda[t1_idx[i]] = exp(log_W[i]+omega);
	}
}
// 
model {
	////////////////////////////////////////////////////////////////////// Bernuli Model
	Z_qst ~ bernoulli_logit(psi_st); ///////////////////////////////////// Standards
	Z_qen_wat ~ bernoulli_logit(psi_un_wat); ///////////////////////////// Water samples
	Z_qen_air ~ bernoulli_logit(psi_un_air); ///////////////////////////// Air samples
	// 
	////////////////////////////////////////////////////////////////////// Continuous Model
	R_qst ~ normal(mu_st,sigma_st); //////////////////////////////////////// Standards
	R_qen_wat ~ normal(mu_en_wat,sigma_en_wat); //////////////////////////// Water samples
	R_qen_air ~ normal(mu_en_air,sigma_en_air); //////////////////////////// Air samples
	//
	////////////////////////////////////////////////////////////////////// Poisson model
	for(i in 1:time){
		N[i] ~ poisson(lambda[i]*E[i]);
	}
	////////////////////////////////////////////////////////////////////// Priors
	//////////////////////////////////////////////////////////////////////// Bernoulli model
	logit_phi ~ normal(logit_phi_mu,logit_phi_sd);
	//////////////////////////////////////////////////////////////////////// Continious model
	beta_0 ~ normal(40, 1);
	beta_1 ~ normal(-3, 1);
	gamma_0 ~ normal(1, 0.1);
	gamma_1 ~ normal(0, 0.1);
  //////////////////////////////////////////////////////////////////////// Count model
	lambda ~ gamma(10,1);
	//////////////////////////////////////////////////////////////////////// Count to Water
	omega ~ normal(1, 1);
	// log_W ~ normal(0,3);
	//////////////////////////////////////////////////////////////////////// Water to Air
	for(i in 1:N_filt){
		eta[i] ~ normal(0,4);
		// eta ~ normal(0,4);
	}
	for(i in 1:N_ft_air){
	epsilon[i] ~ normal(0, exp(tau[air_filt_idx2[i]]));
	// epsilon[i] ~ normal(0, tau[air_filt_idx2[i]]);
	// epsilon ~ normal(0, 1);
	}
	// eta ~ normal(9, 2);
	tau ~ normal(tau_p1,tau_p2);
	log_A ~ normal(0,3);
}

