#' The 'rsstap' package.
#'
#' @description Spline implimentation of Spatial Temporal Aggregated Predictors
#'
#' @docType package
#' @name rsstap-package
#' @aliases rsstap
#' @useDynLib rsstap, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom stats as.formula cov cov2cor gaussian binomial poisson quantile terms
#'
#' @section STAP models:
#' For link function \eqn{g()}, outcome \eqn{\mu_i},  \eqn{i=1,...,n}, STAP Models are 
#' defined by aggregating exposure to nearby built environment features (BEFs) through 
#' either measured space or time.
#' 
#' \deqn{ g(\mu_i) = Z_i^T \beta + \sum f(d) }
#' 
#' 
#' \deqn{ g(\mu_i) = Z_i^T \beta + \sum f(t) }
#' 
#' 
#' \deqn{ g(\mu_i) = Z_i^T \beta + \sum f(d,t) }
#' 
#' Incorporating these models into a standard regression framework is accomplished in this package
#' via a generalized additive approach, using aggregated smooth terms via the \code{\link[mgcv]{jagam}} function.
#' 
#' @section Basic Use:
#' Inspired by the \pkg{mgcv} package, \pkg{rsstap} incorporates STAP terms in the model via the use of a coded
#' term in the standard model formula. For example to fit a model with continuous outcome BMI, covariate sex, and 
#' include distances to FFR as a measure of \emph{spatial}  exposure (not temporal) one would use the following 
#' model formula in (say) the sstap_lm function:
#' 
#' \code{sstap_lm(BMI ~ sex + sap(FFR),benvo = FFR_benvo).}
#' 
#' If aggregating (only) temporal exposure one would use the keyword \code{tap} and if both space and 
#' time are measured, one could use the keyword \code{stap}. The dimension of the basis function expansion
#' can be set familiarly, as a possible option: \code{sap(FFR, k = 5)}. For now the choice of splines is 
#' fixed within this package to be penalized b-splines. This may change in the future.
#'
#' Note that these terms are parsed using regex style functions and so care must be taken not to misspell terms.
#' Any use outside of what is defined here will likely result in erroneous behavior.
#'
#' @section mgcv translation:
#'
#' A table showing the exact rsstap's mgcv counterparts are shown below:
#' \tabular{rr}{
#'   rsstap  \tab mgcv \cr
#'   sap(foo,k=1) \tab s(Distance,k=1,bs='ps')  \cr
#'   sap(foo) \tab s(Distance,bs='ps')  \cr
#'   tap(foo) \tab s(Time,bs='ps')  \cr
#'   stap(foo) \tab t2(Distance,Time,bs='ps') \cr 
#' }
#' 
#' 
#' @references
#' \enumerate{
#'	   \item Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#'     \item  Wood, S.N. (2017) Generalized Additive Models: An Introduction with R (2nd edition).Chapman and Hall/CRC.
#' }
#' @seealso \href{https://apeterson91.github.io/rsstap/articles/Introduction.html}{Introductory Vignette}
NULL
