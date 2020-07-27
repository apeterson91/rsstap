
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

