
#' Small benvo for use in \pkg{rsstap} examples and vignettes.
#'
#' @name binomial_benvo 
#'
#' @format  A benvo with 300 subjects, binomial outcome and nearby simulated "Fast Food Restaurants"(FFRs) and "Healthy Food Stores"(HFS)
#' \describe{
#' \item{\code{subject_data}}{includes Obese and sex binary variables and continuous Income variable}
#' \item{\code{HFS data}}{simulated distances in (0,1) for "nearby HFS"}
#' \item{\code{FFR data}}{simulated distances in (0,1) for "nearby FFRs"}
#' }
#'@seealso the data generation code at <https://github.com/apeterson91/rsstap/tree/master/R/data-raw/binomial.R>
#' 
"binomial_benvo"


#' Small benvo for use in \pkg{rsstap} examples and vignettes.
#'
#' @name poisson_benvo 
#'
#' @format  A benvo with 300 subjects, poisson outcome and nearby simulated "Fast Food Restaurants"(FFRs) and "Healthy Food Stores"(HFS)
#' \describe{
#' \item{\code{subject_data}}{includes Obese count, sex binary variable and continuous Income variable}
#' \item{\code{HFS data}}{simulated distances in (0,1) for "nearby HFS"}
#' \item{\code{FFR data}}{simulated distances in (0,1) for "nearby FFRs"}
#' }
#'@seealso the data generation code at <https://github.com/apeterson91/rsstap/tree/master/R/data-raw/poisson_benvo.R>
#' 
"poisson_benvo"


#' Complex longitudinal benvo for use in \pkg{rsstap} examples and vignettes
#'
#' @name complex_longitudinal
#'
#' @format a benvo with 600 subjects, continuous outcome simulated using subject specific effect and temporal slope and between/within subject exposure effects from nearby "Healthy Food Stores"(HFS)
#' \describe{
#' \item{\code{subject_data}}{Includes pseudo Obesity outcome along with simulated exposure effects}
#' \item{ HFS data  }{ distances and times that subjects spend near HFS }
#' }
#'
#' @seealso The data generation code on \href{github}{https://github.com/apeterson91/rsstap/tree/master/R/data-raw/complex_longitudinal.R}
"complex_longitudinal"


#' Network Built Environment Feature Effects  dataset
#'
#' @name network_benvo
#' 
#' @format a benvo with 500 subjects, a pseudo BMI continuous outcome, pseudo continuous income and binary sex covariate measures with Fast Food Restaurant (FFR)
#' direct exposure effects and "indirect" network effects.
#' 
"network_benvo"
