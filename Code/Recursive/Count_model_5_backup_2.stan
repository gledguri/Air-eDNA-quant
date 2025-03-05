data {
	//Numbers of dimentions
	int N_st_q; // Total number of observation in qPCR standard samples
	int N_en_wat_q; // Total number of observation in qPCR environmental samples
	int N_en_air_q; // Total number of observation in qPCR environmental samples
	// 
		int N_st_qp; // Total number of observation in qPCR standard samples for only detected samples
	int N_en_wat_qp; // Total number of observation in qPCR environmental samples for only detected samples
	int N_en_air_qp; // Total number of observation in qPCR environmental samples for only detected samples
	// 
		int N_filt; // Total number of plates used
	int N_plate; // Total number of plates used
	int N_j_wat; // Number of samples in environmental data
	int N_j_air; // Number of samples in environmental data
	//
		// // Indexes
	// // // Plate index
	array[N_st_qp] int plate_st_idx; // Plate index
	array[N_en_wat_qp] int plate_en_wat_idx; // Plate index for water samples
	array[N_en_air_qp] int plate_en_air_idx; // Plate index for air samples
	// // // Binomial model
	array[N_en_wat_q] int j_qen_wat_idx; // Species and standard index for qPCR environmental samples
	array[N_en_air_q] int j_qen_air_idx; // Species and standard index for qPCR environmental samples
	// // // Continious model
	array[N_en_wat_qp] int j_qen_wat_p_idx; // Species and standard index for qPCR environmental samples
	array[N_en_air_qp] int j_qen_air_p_idx; // Species and standard index for qPCR environmental samples
	// Data
	// // // Binomial model
	array[N_st_q] int Z_qst; // Presence/Absence response of qPCR standard data
	array[N_en_wat_q] int Z_qen_wat; // Presence/Absence response of qPCR environmental data
	array[N_en_air_q] int Z_qen_air; // Presence/Absence response of qPCR environmental data
	// // // Continious model
	array[N_st_qp] real R_qst; // Ct values of qPCR standard data for only detected samples
	array[N_en_wat_qp] real R_qen_wat; // Ct values of qPCR environmental data for only detected samples
	array[N_en_air_qp] real R_qen_air; // Ct values of qPCR environmental data for only detected samples
	// // // 
		array[N_st_q] real S_q; // Known concentration (log10) in qPCR data
	array[N_st_qp] real S_q_p; // Known concentration (log10) in qPCR data for only detected samples
	// 
		int time; 
	array[time] int<lower=0> y;
	array[time] int<lower=0> days;
	array[N_j_wat] int t1_idx;
	array[N_j_wat] int t0_idx;
	array[N_j_air] int air_wat_idx;
	array[N_j_air] int air_filt_idx;
	real tau_p1;
	real tau_p2;
}
// 
	parameters {
		////////////////////////////////////////////////////////////////////// Bernuli parameters
		real alpha_0;
		real alpha_1;
		////////////////////////////////////////////////////////////////////// Continous parameters
		vector[N_plate] eta_0;
		real eta_1;
		real gamma_0;
		real<upper=0> gamma_1;
		////////////////////////////////////////////////////////////////////// Count parameters
		array[time] real<lower=0> lambda;
		////////////////////////////////////////////////////////////////////// Count to Water param
		real beta_0;
		real beta_1;
		////////////////////////////////////////////////////////////////////// Water to Air param
		// array[N_filt] real<lower=0> tau;
		array[N_j_air] real epsilon;
		array[N_filt] real<lower=0,upper=1> eta;
		// Other
		// real<lower=0,upper=1> eta;
		array[N_j_air] real<lower=0> tau;
		// vector[N_j_air] D_air;
	}
// 
	transformed parameters{
		////////////////////////////////////////////////////////////////////// Transf param declare
		// Bernoulli
		vector[N_st_q] theta_st;
		vector[N_en_wat_q] theta_un_wat;
		vector[N_en_air_q] theta_un_air;
		// Continious
		vector[N_st_qp] mu_st;
		vector[N_st_qp] sigma_st;
		vector[N_en_wat_qp] mu_en_wat;
		vector[N_en_wat_qp] sigma_en_wat;
		vector[N_en_air_qp] mu_en_air;
		vector[N_en_air_qp] sigma_en_air;
		// Count to Water
		vector[N_j_wat] D_wat;
		// Water to Air
		vector[N_j_air] D_air;
		// 
			////////////////////////////////////////////////////////////////////// Count to Water
		for (i in 1:N_j_wat){
			// exp(D_wat[i]) = beta_0*lambda[t1_idx[i]]+beta_1*lambda[t0_idx[i]];
			// D_wat[i] = (beta_0*log(lambda[t1_idx[i]])+beta_1*log(lambda[t0_idx[i]]))*(log10(2)/log(2));
			D_wat[i] = (beta_0*log(lambda[t1_idx[i]]))*(log10(2)/log(2));
		}
		////////////////////////////////////////////////////////////////////// Water to Air
		for (i in 1:N_j_air){
			D_air[i] = log(eta[air_filt_idx[i]]*exp(D_wat[air_wat_idx[i]]) + exp(epsilon[i]));
			// D_air[i] = log(eta[air_filt_idx[i]]*exp(D_wat[air_wat_idx[i]]));
			// D_air[i] = (eta[air_filt_idx[i]]*D_wat[air_wat_idx[i]]) + epsilon[i];
		}
		////////////////////////////////////////////////////////////////////// Bernuli Model
		//////////////////////////////////////////////////////////////////////// Standards
		for (i in 1:N_st_q){ 
			theta_st[i] = alpha_0 + (alpha_1 * S_q[i]); 
		}
		//////////////////////////////////////////////////////////////////////// Water samples
		for (i in 1:N_en_wat_q){
			theta_un_wat[i] = alpha_0 + (alpha_1 * D_wat[j_qen_wat_idx[i]]);
		}
		//////////////////////////////////////////////////////////////////////// Air samples
		for (i in 1:N_en_air_q){
			theta_un_air[i] = alpha_0 + (alpha_1 * D_air[j_qen_air_idx[i]]);
		}
		////////////////////////////////////////////////////////////////////// Continuous Model
		//////////////////////////////////////////////////////////////////////// Standards
		for (i in 1:N_st_qp){
			mu_st[i] = eta_0[plate_st_idx[i]] + (eta_1 * S_q_p[i]);
			sigma_st[i] = exp(gamma_0+(gamma_1 * S_q_p[i]));
		}
		//////////////////////////////////////////////////////////////////////// Water samples
		for (i in 1:N_en_wat_qp){
			mu_en_wat[i] = eta_0[plate_en_wat_idx[i]] + (eta_1 * D_wat[j_qen_wat_p_idx[i]]);
			sigma_en_wat[i] = exp(gamma_0+(gamma_1 * D_wat[j_qen_wat_p_idx[i]]));
		}
		//////////////////////////////////////////////////////////////////////// Air samples
		for (i in 1:N_en_air_qp){
			mu_en_air[i] = eta_0[plate_en_air_idx[i]] + (eta_1 * D_air[j_qen_air_p_idx[i]]);
			sigma_en_air[i] = exp(gamma_0+(gamma_1 * D_air[j_qen_air_p_idx[i]]));
		}
	}
// 
	model {
		////////////////////////////////////////////////////////////////////// Bernuli Model
		Z_qst ~ bernoulli_logit(theta_st); ///////////////////////////////////// Standards
		Z_qen_wat ~ bernoulli_logit(theta_un_wat); ///////////////////////////// Water samples
		Z_qen_air ~ bernoulli_logit(theta_un_air); ///////////////////////////// Air samples
		// 
			////////////////////////////////////////////////////////////////////// Continuous Model
		R_qst ~ normal(mu_st,sigma_st); //////////////////////////////////////// Standards
		R_qen_wat ~ normal(mu_en_wat,sigma_en_wat); //////////////////////////// Water samples
		R_qen_air ~ normal(mu_en_air,sigma_en_air); //////////////////////////// Air samples
		// 
			////////////////////////////////////////////////////////////////////// Poisson model
		for(i in 1:time){
			y[i] ~ poisson(lambda[i]*days[i]);
		}
		////////////////////////////////////////////////////////////////////// Priors
		//////////////////////////////////////////////////////////////////////// Bernoulli model
		alpha_0 ~ normal(0, 1);
		alpha_1 ~ normal(0, 1);
		//////////////////////////////////////////////////////////////////////// Continious model
		eta_0 ~ normal(40, 1);
		eta_1 ~ normal(-3, 1);
		gamma_0 ~ normal(1, 0.1);
		gamma_1 ~ normal(0, 0.1);
		//////////////////////////////////////////////////////////////////////// Count model
		lambda ~ gamma(10,1);
		//////////////////////////////////////////////////////////////////////// Count to Water
		beta_0 ~ normal(0, 2);
		beta_1 ~ normal(0, 2);
		// D_wat ~ normal(0,3);
		//////////////////////////////////////////////////////////////////////// Water to Air
		// for(i in 1:N_j_air){
			// 	epsilon[i] ~ normal(0, tau[air_filt_idx[i]]);
			// }
		tau ~ gamma(tau_p1,tau_p2);
		eta ~ beta(1, 1);
		// Other
		epsilon ~ normal(0, tau);
		// D_air ~ normal(0,3);
	}

