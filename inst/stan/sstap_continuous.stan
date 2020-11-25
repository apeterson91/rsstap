// GLM for a Gaussian outcome with smooth terms
functions {

#include functions/common_functions.stan
  
}
data{ 
	int<lower=1> N;
	int<lower=1> ncol_Z;
	int<lower=1> ncol_smooth;
	int<lower=1> K_smooth;
	int<lower=1> num_stap;
	int<lower=1> stap_lengths[num_stap];
	int<lower=1> stap_penalties[num_stap];
	int<lower=1> num_stap_penalties;
	int<lower=0> stap_pen_map[num_stap,max(stap_penalties)];
	int<lower=1> pen_ix[num_stap_penalties,2];
	int<lower=1> beta_ix[num_stap,2];
	int<lower=1> P;
	vector[N] y;
	matrix[N,P] Q;
	matrix[P,P] R_inv;
	matrix[ncol_smooth,K_smooth] S;

	//data for glmer
#include data/glmer_stuff.stan
#include data/glmer_stuff2.stan

}
transformed data{
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);
#include tdata/tdata_glmer.stan
}
parameters{
	vector[P] beta_tilde;
	vector<lower=0>[num_stap_penalties] tau;
	real<lower=0> sigma;
	//merModels
#include parameters/parameters_glmer.stan
}
transformed parameters{
	vector[P] beta = R_inv * beta_tilde;
	vector[ncol_smooth] sstap_beta = beta[(ncol_Z+1):P];
	vector[N] eta;
	//merModels
	// construct L, b
#include tparameters/tparameters_glmer.stan

	//construct eta
	eta = Q * beta_tilde;

	if(t>0){
		if (special_case) for (i in 1:t) eta += b[V[i]];
		else eta += csr_matrix_times_vector(N, q, w, v, u, b);
	}


}
model{

	sigma ~ cauchy(0,5);
	tau ~ exponential(1);
	y ~ normal(eta,sigma);


	// add penalty
	for(i in 1:num_stap){
		matrix[stap_lengths[i],stap_lengths[i]] K = rep_matrix(0,stap_lengths[i],stap_lengths[i]);
		for(j in 1:stap_penalties[i]){
			int pen_mapper = stap_pen_map[i,j];
			K = K + S[ beta_ix[i,1]:beta_ix[i,2], pen_ix[pen_mapper,1]:pen_ix[pen_mapper,2]] * tau[pen_mapper];
		}
		target += multi_normal_prec_lpdf(sstap_beta[ beta_ix[i,1]:beta_ix[i,2]  ]|  rep_vector(0,stap_lengths[i]),K );
	}

  if (t > 0) {
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau_b, 
                          regularization, delta, shape, t, p);
  }


}
generated quantities {
	vector[ncol_Z] delta_coef = beta[1:ncol_Z];
}
