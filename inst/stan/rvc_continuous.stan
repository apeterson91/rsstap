//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
functions{
     vector pw_gaussian(vector y, vector eta, real sigma) {
    return -0.5 * log(6.283185307179586232 * sigma) -
            0.5 * square((y - eta) / sigma);
  }
}
data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> L;
  int<lower=0> M;
  int<lower=0> K;
  int<lower=0> w_num;
  int<lower=0> v_num;
  int<lower=0> v[v_num];
  vector[w_num] w;
  int<lower=0> u[N+1];
  vector[N] y;
  matrix[M,L] X_smooth;
  matrix[N,P] X_fixef;
  matrix[L,L] S1;
  matrix[L,L] S2;
  matrix[M,K] D;
	// data for weights
#include data/weights.stan

}
parameters {
  vector[P] beta_fixef;
  vector[L] beta_smooth;
  vector[K] alpha;
  real<lower=0> sigma;
  real<lower=0> tau1;
  real<lower=0> tau2;
}

model {
  vector[M] eta;
  vector[N] eta_;
  alpha ~ normal(0,1);
  tau1 ~ gamma(1,1);
  tau2 ~ gamma(1,1);
  sigma ~ cauchy(0,5);
  beta_smooth ~ multi_normal_prec(rep_vector(0,L),S1 * tau1 + S2 * tau2);
  eta = exp(D*alpha) .* (X_smooth * beta_smooth);
  eta_ = csr_matrix_times_vector(N,M,w,v,u,eta) + X_fixef * beta_fixef;
  if(has_weights)
	  target += dot_product(weights,pw_gaussian(y,eta_,sigma));
  else
	  y ~ normal(eta_,sigma);
}
