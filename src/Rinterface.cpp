// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
//


#include "BData.hpp"
#include "Tree.hpp"
#include "NUTS.hpp"


//' Fits a bayesian built environment network regression model using the No-U-Turn Sampler
//'
//' @param y vector of outcomes
//' @param Z matrix of subject covariates
//' @param DD sparse matrix of subject distances in basis format
//' @param rpc array of regression parameter settings
// [[Rcpp::export]]
Rcpp::List bnet_lm_fit(const Eigen::VectorXd &y,
                       const Eigen::MatrixXd &Z,
                       const SEXP &DD,
                       const Eigen::ArrayXi &rpc,
                       const Eigen::ArrayXd &prior_means,
                       const Eigen::ArrayXd &prior_scales,
                       const int &iter_max,
                       const int &warm_up,
                       const int &seed) {

	auto start = std::chrono::high_resolution_clock::now();

	const Eigen::MappedSparseMatrix<double> D(Rcpp::as<Eigen::MappedSparseMatrix<double> >(DD));

	const Bdata(y,Z,D,rpc,prior_means,prior_scales);

	std::mt19937 rng;
	rng = std::mt19937(seed);

	for(int iter_ix = 1; iter_ix ++ ; iter_ix < iter_max){


	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
	double sampling_time = duration.count();

    return(Rcpp::List::create(Rcpp::Named("sampling_time")=0);
}
