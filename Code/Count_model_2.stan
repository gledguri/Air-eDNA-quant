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
  /**
  * Compute zero-sum constrained parameter
  *
  * @param delta_raw A declared parameter in Stan
  * @param tau_raw A declared parameter in Stan
  
  * @param N_samp Number of samples
  * @param N_samp_bio_rep Number of samples and bio_replicate
  * @param N_samp_w_bio_rep Number of samples with bio_replicate
  
  * @param bio_rep_idx Vector of the total number of bio_replicate per sample
  * @param tau_bio_rep_idx Index vector for having idiosyncratic tau parameter
  
  * @return Vector[N_samp_bio_rep] with delta (sum-zero constrained) parameters
  */
  vector sum_zero_constrain(vector delta_raw, 
                       vector tau_raw, 
                       int N_samp,
                       int N_samp_bio_rep,
                       array[] int bio_rep_idx, 
                       array[] int tau_bio_rep_idx) {
    
    vector[N_samp_bio_rep] delta;
    int count_tot = 0;
    int count_par = 0;
    real delta_sum;
    
    for (j in 1:N_samp) {
      delta_sum = 0;
      
      for (k in 1:bio_rep_idx[j]) {
        count_tot = count_tot + 1;
        
        if (k < bio_rep_idx[j]) {
          count_par = count_par + 1;
          delta[count_tot] = delta_raw[count_par] * tau_raw[tau_bio_rep_idx[j]];
          delta_sum = delta_sum + delta[count_tot];
        } else if (bio_rep_idx[j] == 1) {
          delta[count_tot] = 0;
        } else {
          delta[count_tot] = -delta_sum;
        }
      }
    }
    
    return delta;
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
// 	// // // Binomial mode
	array[N_st_q] int Z_qst; // Presence/Absence response of qPCR standard data
	array[N_en_wat_q] int Z_qen_wat; // Presence/Absence response of qPCR environmental data
	array[N_en_air_q] int Z_qen_air; // Presence/Absence response of qPCR environmental data
// 	// // // Continious model
	vector[N_st_qp] R_qst; // Ct values of qPCR standard data for only detected samples
	vector[N_en_wat_qp] R_qen_wat; // Ct values of qPCR environmental data for only detected samples
	vector[N_en_air_qp] R_qen_air; // Ct values of qPCR environmental data for only detected samples
// 	// // // 
	vector[N_st_q] S_q; // Known concentration (log10) in qPCR data
	vector[N_st_qp] S_q_p; // Known concentration (log10) in qPCR data for only detected samples
// 	// 
	int time;
	array[time] int<lower=0> N;
	// array[time] int<lower=0> E;
	vector[time] E;
	// 
	// array[N_j_wat] int t_lambda_idx;
	array[time] int l_i;
	// 
	array[N_j_air] int a_j;
	array[N_j_air] int a_i;
	array[N_j_air] int a_ij;
	array[N_j_air] int a_ijb;
	array[N_j_air] int bio_rep;
	array[N_j_wat] int W_i_idx;
	// 
	real tau_p1;
	real tau_p2;
	real logit_phi_mu;
	real logit_phi_sd;
	int N_samp_bio_rep;
  int N_samp_w_bio_rep; 
  int N_samp; 
  array[N_samp] int bio_rep_idx; 
  array[N_samp] int tau_bio_rep_idx; 
}
parameters {
	////////////////////////////////////////////////////////////////////// Bernuli parameters
	real logit_phi;
	////////////////////////////////////////////////////////////////////// Continous parameters
  vector[N_plate] beta_0;
	real beta_1;
	real gamma_0;
	real<upper=0> gamma_1;
	////////////////////////////////////////////////////////////////////// Count to Water param
	vector[N_filt] eta;
	vector[N_ft_air] epsilon_raw;
	////////////////////////////////////////////////////////////////////// Water to X
	real omega;
	////////////////////////////////////////////////////////////////////// X parameter
	vector<lower=0>[N_j_wat] X_STATE;
	vector<lower=0>[N_filt] tau;
	vector[N_samp_w_bio_rep] delta_raw; 
	vector<lower=0>[N_filt] tau_raw;
	// real theta;
	// vector[N_j_wat] log_W;
	// vector[N_j_wat] kappa_raw;
	// real<lower=0> kappa_sd;
  // real<lower=0> tau_raw; 
  // vector[N_j_wat] kappa;
}
// 
transformed parameters{
  vector[N_samp_bio_rep] delta = sum_zero_constrain(delta_raw,tau_raw,N_samp,N_samp_bio_rep,bio_rep_idx,tau_bio_rep_idx);
	
	////////////////////////////////////////////////////////////////////// Count to Water
	// vector[N_j_wat] kappa = kappa_raw*kappa_sd;
	vector[N_j_wat] log_W = log(X_STATE)-omega;	
	//////////////////////////////////////////////////////////////////////// Air to Water
	vector[N_ft_air] epsilon = epsilon_raw .* tau[tau_bio_rep_idx];
	vector[N_j_air] log_A;
	log_A[a_ijb] = eta[a_j] + (log_W[a_i]) + epsilon[a_ij] + delta[a_ijb];
	//////////////////////////////////////////////////////////////////////// Standards
	vector[N_st_q] p_tmp_st = exp(S_q) * -inv_logit(logit_phi);
  vector[N_st_q] psi_st = log1m_exp(p_tmp_st) - p_tmp_st;
	vector[N_st_qp] mu_st = beta_0[plate_st_idx] + (beta_1 * S_q_p);
	vector[N_st_qp] sigma_st = exp(gamma_0+(gamma_1 * S_q_p));
	//////////////////////////////////////////////////////////////////////// Water samples
  vector[N_en_wat_q] p_tmp_un_wat = exp(log_W[j_qen_wat_idx]) * -inv_logit(logit_phi);
  vector[N_en_wat_q] psi_un_wat = log1m_exp(p_tmp_un_wat) - p_tmp_un_wat;
	vector[N_en_wat_qp] mu_en_wat = beta_0[plate_en_wat_idx] + (beta_1 * log_W[j_qen_wat_p_idx]);
	vector[N_en_wat_qp] sigma_en_wat = exp(gamma_0+(gamma_1 * log_W[j_qen_wat_p_idx]));
	//////////////////////////////////////////////////////////////////////// Air samples
	vector[N_en_air_q] p_tmp_un_air =  exp(log_A[j_qen_air_idx]) * -inv_logit(logit_phi);
  vector[N_en_air_q] psi_un_air = log1m_exp(p_tmp_un_air) - p_tmp_un_air;
	vector[N_en_air_qp] mu_en_air = beta_0[plate_en_air_idx] + (beta_1 * log_A[j_qen_air_p_idx]);
	vector[N_en_air_qp] sigma_en_air = exp(gamma_0+(gamma_1 * log_A[j_qen_air_p_idx]));
	vector[time] lambda = X_STATE[l_i].*E;
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
	N ~ neg_binomial_2(lambda,20);
	// N ~ poisson(lambda);
	////////////////////////////////////////////////////////////////////// Priors
	//////////////////////////////////////////////////////////////////////// Bernoulli model
	logit_phi ~ normal(logit_phi_mu,logit_phi_sd);
	//////////////////////////////////////////////////////////////////////// Continious model
	beta_0 ~ normal(40, 1);
	beta_1 ~ normal(-3, 1);
	gamma_0 ~ normal(1, 0.1);
	gamma_1 ~ normal(0, 0.1);
	// theta ~ normal(1, 0.1);
  //////////////////////////////////////////////////////////////////////// Count model
	X_STATE ~ normal(100,100);
	//////////////////////////////////////////////////////////////////////// Count to Water
	omega ~ normal(-6, 1);
	// log_W ~ normal(0,3);
	//////////////////////////////////////////////////////////////////////// Water to Air
	log_A ~ normal(0,3);
	eta ~ normal(-2,4);
	epsilon_raw ~ std_normal();
	tau ~ std_normal();
	delta_raw ~ std_normal(); 
  tau_raw ~ std_normal();
}
