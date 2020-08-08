

#' Print method for sstapreg objects
#'
#' The \code{print} method for stanreg objects displays a compact summary of the
#' fitted model. 
#' @export
#' @method print sstapreg
#' @param x sstapreg object
#' @param digits number of digits to round to
#' @param ... ignored
#'
print.sstapreg <- function(x, digits = 1, ...) {

    cat("\n family:      ", family_plus_link(x))
    cat("\n formula:     ", formula_string(formula(x)))
    cat("\n observations:", nobs(x))
	cat("\n fixed effects: ", length(coef(x)))
	if(any(x$stap_components == "Distance"))
		cat("\n spatial predictors: ", sum(x$stap_components == "Distance"))
	if(any(x$stap_components == "Time"))
		cat("\n temporal predictors: ", sum(x$stap_components == "Time"))
	if(any(x$stap_components == "Distance-Time"))
		cat("\n spatial-temporal predictors: ", sum(x$stap_components == "Distance-Time"))
    cat("\n------\n")

	mat <- cbind(coef(x),se(x))
	colnames(mat) <- c("Median","MAD")

	.printfr(mat,digits)

    mat <- as.matrix(x$stapfit) # don't used as.matrix.stanreg method b/c want access to mean_PPD

	smooth_nms <- grep("^smooth_precision\\[", colnames(mat), value = TRUE)
	smooth_mat <- mat[,smooth_nms,drop=F]
	smooth_ests <- .median_and_madsd(smooth_mat)

      cat("\nSmoothing terms:\n")
	.printfr(smooth_ests,digits)
    if (is.mer(x)) {
      cat("\nError terms:\n")
      foo <- VarCorr(x)
	  print(foo)
      cat("Num. levels:", 
          paste(names(ngrps(x)), unname(ngrps(x)), collapse = ", "), "\n")
    }



}

# Family name with link in parenthesis 
# @param x sstapreg object
family_plus_link <- function(x) {
  fam <- family(x)
  fam <- paste0(fam$family, " [", fam$link, "]")
  return(fam)
}

# @param formula formula object
formula_string <- function(formula, break_and_indent = TRUE) {
  coll <- if (break_and_indent) "--MARK--" else " "
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = coll))
  if (!break_and_indent)
    return(char)
  gsub("--MARK--", "\n\t  ", char, fixed = TRUE)
}

.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}
.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}
