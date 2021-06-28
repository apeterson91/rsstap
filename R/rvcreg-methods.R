#' Methods for rvcreg objects
#' 
#' The methods documented on this page are actually some of the least important 
#' methods defined for \link{rvcreg} objects. The most 
#' important methods are documented separately, each with its own page. Links to
#' those pages are provided in the \strong{See Also} section, below.
#' This documentation is inspired by the \code{\link[rstanarm]{rstanarm}} package
#' 
#' @name rvcreg-methods
#' 
#' @param object rvcreg object
#' @param ... Ignored
#' 
#' @importFrom stats coef nobs 
#' @details The methods documented on this page are similar to the methods 
#'   defined for objects of class 'lm', 'glm', 'glmer', etc. However there are a
#'   few key differences:
#'   
#' \describe{
#' \item{\code{coef}}{
#' Medians are used for fixed effect point estimates. See the \emph{Point estimates} section
#' in \code{\link{print.sstapreg}} for more details.
#' }
#' \item{\code{se}}{
#' The \code{se} function returns fixed effects standard errors based on 
#' \code{\link{mad}}. See the \emph{Uncertainty estimates} section in
#' \code{\link{print.sstapreg}} for more details.
#' }
#' \item{\code{confint}}{
#' \code{confint} will throw an error because the
#' \code{\link{posterior_interval}} function should be used to compute Bayesian 
#' uncertainty intervals.
#' }}
#' 
#' 
#' 
#' 
NULL

#' @rdname rvcreg-methods
#' @export
nobs.rvcreg <- function(x)
	return(length(x$spec$mf$y))

#' @rdname rvcreg-methods
#' @export
coef.rvcreg <- function(x)
	return(x$coefficients)

#' @rdname rvcreg-methods
#' @export
ses.rvcreg <- function(x)
	return(x$ses)

#' @rdname rvcreg-methods
#' @export
confint.rvcreg <- function(object, ...) {
    stop("Please use posterior_interval() to obtain ", 
         "Bayesian interval estimates.", 
         call. = FALSE)
}

#' family method for rvcreg objects
#' 
#' @keywords internal
#' @export 
#' @param object,... See \code{\link[stats]{family}}
family.rvcreg <- function(object,...) object$family

#' @rdname rvcreg-methods
#' @export
se.rvcreg <- function(object, ...) {
  object$ses
}
