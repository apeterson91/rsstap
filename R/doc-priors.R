#' Priors for `rsstap` models
#'
#' @name priors
#'
#' @section GLMs: 
#' Using  one spatial aggregated predictor as an example, \code{stap_(g)lm}  models have the following form.
#' 
#'
#' \deqn{ g(\mu_i) = Z_i^T \delta + \sum_d \sum_l \beta_l\phi_l(d) }
#'
#'
#'
#' Currently the priors in \pkg{rsstap} are fixed and are always of the following form:
#'
#'  \deqn{p(\delta) \propto  1}
#'
#'  \deqn{\sigma \sim C^+(0,5)}
#' 
#'  \deqn{\beta \sim MVN_L(0,\sum_k S_k \tau_k)}
#'
#'  \deqn{\tau_k \sim \text{Exponential}(1)}
#'
#' Where \eqn{S_k} are generated from the \code{\link[mgcv]{jagam}} function and sum to form a complete precision matrix with different \eqn{\tau} penalties along the diagonal.
#'
#' @section GLMERs:
#'
#' Using only one spatial aggregated predictor as an example,  \code{stap_(g)lmer}  models have the following form:
#'
#' \deqn{ g(\mu_{ij}) = Z_{ij}^T \delta + \sum_d \sum_l \beta_l\phi_l(d) + W_{ij}^Tb_i }
#'
#' Where
#' 
#' \deqn{b_i \sim N(0,\Sigma)}
#'
#' priors for \eqn{\delta},\eqn{\beta},\eqn{\sigma},\eqn{\tau_k} are the same as before, but now  
#' \eqn{\Sigma} is decomposed as described \href{https://mc-stan.org/rstanarm/articles/glmer.html#priors-on-covariance-matrices-1}{here}.
#' 
NULL
