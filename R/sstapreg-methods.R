#' Methods for sstapreg objects
#' 
#' The methods documented on this page are actually some of the least important 
#' methods defined for \link[=sstapreg-objects]{sstapreg} objects. The most 
#' important methods are documented separately, each with its own page. Links to
#' those pages are provided in the \strong{See Also} section, below.
#' 
#' @name sstapreg-methods
#' @aliases VarCorr fixef sigma
#' 
#' @param object sstapreg object
#' @param ... Ignored
#' 
#' @details The methods documented on this page are similar to the methods 
#'   defined for objects of class 'lm', 'glm', 'glmer', etc. However there are a
#'   few key differences:
#'   
#' \describe{
#' \item{\code{coef}}{
#' Medians are used for point estimates. See the \emph{Point estimates} section
#' in \code{\link{print.sstapreg}} for more details.
#' }
#' \item{\code{se}}{
#' The \code{se} function returns standard errors based on 
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
 object$Z_coefs
}

#' @rdname sstapreg-methods
#' @export
#'
confint.sstapreg <- function(object, ...) {
    stop("Please use posterior_interval() to obtain", 
         "Bayesian interval estimates.", 
         call. = FALSE)
}

#' @rdname sstapreg-methods
#' @export
fitted.sstapreg <- function(object, ...)  
    return(object$fitted.values)

#' @rdname sstapreg-methods
#' @export 
nobs.sstapreg <- function(object, ...) {
  nrow(model.frame(object))
}

#' Extract standard errors
#' 
#' Generic function for extracting standard errors from fitted models.
#' 
#' @export
#' @keywords internal
#' @return Standard errors of model parameters.
#' @seealso \code{\link{se.sstapreg}}
#' 
se <- function(object, ...) UseMethod("se")

#' @rdname sstapreg-methods
#' @export
se.sstapreg <- function(object, ...) {
  object$ses
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


