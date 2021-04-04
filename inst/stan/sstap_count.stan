// GLM for a Poisson outcome with smooth terms
functions {
#include functions/common_functions.stan
#
  /** 
  * Pointwise (pw) log-likelihood vector for the Poisson distribution
  *
  * @param y The integer array corresponding to the outcome variable.
  * @param eta The vector of linear predictors
  * @param link An integer indicating the link function
  * @return A vector
  */
  vector pw_pois(int[] y, vector eta, int link) {
    int N = rows(eta);
    vector[N] ll;
    if (link == 1)  // log
      for (n in 1:N) ll[n] = poisson_log_lpmf(y[n] | eta[n]);
    else reject("Invalid link");
    return ll;
  }

}
data{ 
	int<lower=1> N;
	int<lower=0> y[N];
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
	matrix[N,P] Q;
	matrix[P,P] R_inv;
	matrix[ncol_smooth,K_smooth] S;
	// data for weights
#include data/weights.stan
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
#include parameters/parameters_glmer.stan
}
transformed parameters{
	vector[P] beta = R_inv * beta_tilde;
	vector[ncol_smooth] sstap_beta = beta[(ncol_Z+1):P];
	vector[N] eta;
	//merModels
	// construct L, b
#include tparameters/tparameters_glmer_noaux.stan
	if(t>0){
		if (special_case) for (i in 1:t) eta += b[V[i]];
		else eta += csr_matrix_times_vector(N, q, w, v, u, b);
	}

	eta = Q * beta_tilde;
	eta = exp(eta);

}
model{
	tau ~ exponential(1);
	if(has_weights)
		target += dot_product(weights,pw_pois(y,eta,1)); // Currently only allows log link
	else
		y ~ poisson(eta);


	for(i in 1:num_stap){
		matrix[stap_lengths[i],stap_lengths[i]] K = rep_matrix(0,stap_lengths[i],stap_lengths[i]);
		for(j in 1:stap_penalties[i]){
			int pen_mapper = stap_pen_map[i,j];
			K = K + S[ beta_ix[i,1]:beta_ix[i,2], pen_ix[pen_mapper,1]:pen_ix[pen_mapper,2]] * tau[pen_mapper];
		}
		target += multi_normal_prec_lpdf(sstap_beta[ beta_ix[i,1]:beta_ix[i,2]  ]|  rep_vector(0,stap_lengths[i]),K );
	}


}
generated quantities {
	vector[ncol_Z] delta_coef = beta[1:ncol_Z];
}
