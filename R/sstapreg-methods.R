#' Methods for sstapreg objects
#' 
#' The methods documented on this page are actually some of the least important 
#' methods defined for \link[=sstapreg]{sstapreg} objects. The most 
#' important methods are documented separately, each with its own page. Links to
#' those pages are provided in the \strong{See Also} section, below.
#' 
#' @name sstapreg-methods
#' 
#' @param object sstapreg object
#' @param ... Ignored
#' 
#' @importFrom stats coef nobs formula family
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

#' @rdname sstapreg-methods
#' @export
coef.sstapreg <- function(object, ...) {
 object$coefficients
}

#' @rdname sstapreg-methods
#' @export
#'
confint.sstapreg <- function(object, ...) {
    stop("Please use posterior_interval() to obtain ", 
         "Bayesian interval estimates.", 
         call. = FALSE)
}

#' @rdname sstapreg-methods
#' @export
fitted.sstapreg <- function(object, ...)  
    return(object$fitted.values)

#' Extract standard errors
#' 
#' Generic function for extracting standard errors from fitted models.
#' 
#' @export
#' @keywords internal
#' @param object a Fitted sstapreg model object
#' @param ... argument to method
#' @return Standard errors of model parameters.
#' @seealso \code{\link{se.sstapreg}}
#' 
se <- function(object, ...) UseMethod("se")

#' @rdname sstapreg-methods
#' @export
se.sstapreg <- function(object, ...) {
  object$ses
}


#' @rdname sstapreg-methods
#' @export
nobs.sstapreg <- function(object, ...){
	length(object$fitted.values)
}


#'
#' @rdname sstapreg-methods
#'
#' @export 
#' @param correlation For \code{vcov}, if \code{FALSE} (the default) the
#'   covariance matrix is returned. If \code{TRUE}, the correlation matrix is
#'   returned instead.
#'
vcov.sstapreg <- function(object, correlation = FALSE, ...) {
  out <- object$covmat
  if (!correlation) return(out)
  cov2cor(out)
}

#' @rdname sstapreg-methods
#' @export
#' @export nsamples
#' @importFrom rstantools nsamples
nsamples.sstapreg <- function(object, ...) {
  posterior_sample_size(object)
}

# Exported but doc kept internal ----------------------------------------------

#' family method for sstapreg objects
#' 
#' @keywords internal
#' @export 
#' @param object,... See \code{\link[stats]{family}}
family.sstapreg <- function(object,...) object$family



#' Formula method for sstapreg objects
#'
#' @keywords internal
#' @export
#' @param x a sstapreg object
#' @param ... ignored currently
formula.sstapreg <- function(x,...){
	x$formula
}


