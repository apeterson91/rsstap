// GLM for a Poisson outcome with smooth terms
functions {
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
}
parameters{
	vector[P] beta_tilde;
	vector<lower=0>[num_stap_penalties] tau;
}
transformed parameters{
	vector[P] beta = R_inv * beta_tilde;
	vector[ncol_smooth] sstap_beta = beta[(ncol_Z+1):P];
	vector[N] eta = exp(Q * beta_tilde);
}
model{
	tau ~ exponential(1);
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
	int<lower=0> yhat[N];
	vector[ncol_Z] delta = beta[1:ncol_Z];
	yhat = poisson_rng(eta);
}
