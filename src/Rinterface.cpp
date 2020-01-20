// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <random>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]


void print_progress(const int &iter_ix, const int &warm_up, const int &iter_max, const int &chain){

  if(iter_max > 20){
      if((iter_ix) % (int)round(.1 * iter_max) == 0 || iter_ix == 1 || iter_ix == (warm_up + 1) ){
          int progress = (int)round(iter_ix * 100 / iter_max);
          std::string str = (iter_ix) <= (warm_up) ? "\t [Warmup]" : "\t [Sampling]";
          Rcpp::Rcout << "[Chain " << chain << "] Beginning of iteration: " << (iter_ix) << " / " << iter_max << " (" << progress << "%)" << str  << std::endl;
      }
  }
  else{
          int progress = (int)round(iter_ix * 100 / iter_max);
          std::string str = (iter_ix) <= (warm_up) ? "\t [Warmup]" : "\t [Sampling]";
          Rcpp::Rcout << "[Chain " << chain << "] Beginning of iteration: " << (iter_ix) << " / " << iter_max << " (" << progress << "%)" << str  << std::endl;
  }

}


double initialize_scalar(std::mt19937 &rng){
    
    std::uniform_real_distribution<double> runif(0,1);
    double out = runif(rng);
    
    return(out);
}

Eigen::VectorXd initialize_vec(const int n, std::mt19937 &rng ){
    
    Eigen::VectorXd out(n);
    
    for(int i = 0; i < n ; i ++)
        out(i) = initialize_scalar(rng);
    
    return(out);
}

#include "BData.hpp"
#include "BbNet.hpp"
#include "Tree.hpp"
#include "NUTS.hpp"

//' Bayesian Built Environment Network Model
//' @param y vector of outcomes
//' @param Z matrix of subject covariates
//' @param DD sparse matrix of subject distances in basis format
//' @param rpc array of regression parameter settings
//' @param start_stop array of start-stop beta par indeces for D matrix
//' @param prior_means real valued vector for holding prior mean information
//'  prior_means[0] = prior mean for intercept
//'  prior_means[1] = prior mean for delta coefs ... and so on
//' @param prior_scales real valued vector for holding prior scale information
//'  prior_scales[0] = prior scale for intercept
//'  prior_scales[1] = prior scale for delta coefs ... and so on
//' @param iter_max maximum number of iterations
//' @param max_treedepth max treedepth for NUTS
//' @param warm_up number of iterations to use for epsilon adaptation
//' @param seed random number generator seed
//
// [[Rcpp::export]]
Rcpp::List bbnet_lm_fit(const Eigen::VectorXd &y,
                       const Eigen::MatrixXd &Z,
                       const SEXP DD,
                       const Eigen::ArrayXi &rpc,
					   const Eigen::ArrayXi &start_stop,
                       const Eigen::ArrayXd &prior_means,
                       const Eigen::ArrayXd &prior_scales,
                       const int &iter_max,
					   const int &max_treedepth,
                       const int &warm_up,
                       const int &seed) {

	auto start = std::chrono::high_resolution_clock::now();
	const int chain = 1;

	const Eigen::MappedSparseMatrix<double> D(Rcpp::as<Eigen::MappedSparseMatrix<double> >(DD));

	const BData data(y,Z,D,rpc,start_stop,prior_means,prior_scales);

	std::mt19937 rng;
	rng = std::mt19937(seed);

	for(int iter_ix = 1; iter_ix < iter_max; iter_ix ++ ){
		print_progress(iter_ix, warm_up, iter_max, chain);


	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
	double sampling_time = duration.count();

    return(Rcpp::List::create(Rcpp::Named("sampling_time")=0));

}
