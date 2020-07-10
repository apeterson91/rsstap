

#' Print method for sstapreg objects
#'
#' The \code{print} method for stanreg objects displays a compact summary of the
#' fitted model. See the \strong{Details} section below for descriptions of the
#' different components of the printed output. 
#' @export
#' @method print sstapreg
#' @param x sstapreg object
#' @param digits number of digits to round to
#'
print.sstapreg <- function(x, digits = 1, ...) {

	cat("\n call:", x$call)
}
