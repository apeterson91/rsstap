#' Print method for sstapreg objects
#'
#' The \code{print} method for sstapreg objects displays a compact summary of the
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

    mat <- as.matrix(x$stapfit) # don't used as.matrix.stapreg method b/c want access to mean_PPD

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


#' Summary method for sstapreg objects
#' 
#' Summaries of parameter estimates and MCMC convergence diagnostics 
#' (Monte Carlo error, effective sample size, Rhat).
#'
#' @export
#' @method summary sstapreg 
#'
#' @param object sstapreg object
#' @param probs an optional numeric vector of probabilities pased to \code{\link[stats]{quantile}}
#' @param ... ignored
#' @param digits Number of digits to use for formatting numbers when printing. 
#'   When calling \code{summary}, the value of digits is stored as the 
#'   \code{"print.digits"} attribute of the returned object.
#' @return The \code{summary} method returns an object of class 
#'   \code{"summary.sstapreg"}  which is a matrix of 
#'   summary statistics and diagnostics, with attributes storing information for use by the
#'   \code{print} method. 
#'
#'
#' @importMethodsFrom rstan summary
summary.sstapreg <- function(object,probs=c(0.05,0.5,0.95),...,digits=1){

	args <- list(object = object$stapfit, probs=probs) 
	out <- do.call("summary",args)$summary
	nms <- c(names(coef(object)),
			 grep(pattern="smooth",x=rownames(out),value=T),
			 paste0("b[",make_b_nms(object$glmod$reTrms), "]"),
			 grep(pattern="log-posterior",x=rownames(out),value=T))
	out <- out[nms,]

    stats <- colnames(out)
    if ("n_eff" %in% stats) {
      out[, "n_eff"] <- round(out[, "n_eff"])
    }
    if ("se_mean" %in% stats) {# So people don't confuse se_mean and sd
      colnames(out)[stats %in% "se_mean"] <- "mcse"
    }

	structure(
	out,
	call = object$call,
	family = family_plus_link(object),
	formula = formula(object),
    posterior_sample_size = posterior_sample_size(object),
    nobs = nobs(object),
	ngrps = if(is.mer(object)) ngrps(object) else NULL,
	print.digits = digits,
	class = "summary.sstapreg"
	)
}


#' @rdname summary.sstapreg
#' @export
#' @method print summary.sstapreg
#'
#' @param x An object of class \code{"summary.sstapreg"}.
print.summary.sstapreg <-
  function(x, digits = max(1, attr(x, "print.digits")),
           ...) {

    atts <- attributes(x)
    cat("\nModel Info:")
    cat("\n family:      ", atts$family)
    cat("\n formula:     ", formula_string(atts$formula))
	cat("\n sample:      ", atts$posterior_sample_size, 
	  "(posterior sample size)")
    
    cat("\n observations:", atts$nobs)
    if (!is.null(atts$ngrps)) {
      cat("\n groups:      ", paste0(names(atts$ngrps), " (", 
                                     unname(atts$ngrps), ")", 
                                     collapse = ", "))
    }
	cat("\n")
    
	hat <- "Rhat"
	str_diag <- "MCMC diagnostics"
	str1 <- "and Rhat is the potential scale reduction factor on split chains"
	str2 <- " (at convergence Rhat=1).\n"
    sel <- which(colnames(x) %in% c("mcse", "n_eff", hat))
    has_mc_diagnostic <- length(sel) > 0
    if (has_mc_diagnostic) {
      xtemp <- x[, -sel, drop = FALSE]
      colnames(xtemp) <- paste(" ", colnames(xtemp))
    } else {
      xtemp <- x
    }
    
    # print table of parameter stats
    .printfr(xtemp, digits)
    
    if (has_mc_diagnostic) {
      cat("\n", str_diag, "\n", sep = '')
      mcse_hat <- format(round(x[, c("mcse", hat), drop = FALSE], digits), 
                          nsmall = digits)
      n_eff <- format(x[, "n_eff", drop = FALSE], drop0trailing = TRUE)
      print(cbind(mcse_hat, n_eff), quote = FALSE)
      cat("\nFor each parameter, mcse is Monte Carlo standard error, ", 
          "n_eff is a crude measure of effective sample size, ", 
          str1, 
          str2, sep = '')
    }
    
    invisible(x)
  }

#' Print method for rvcreg objects
#' @export
#' @method print rvcreg 
#' @param x sstapreg object
#' @param digits number of digits to round to
#' @param ... ignored
print.rvcreg <- function(x, digits = 2, ...){

    cat("\n family:      ", family_plus_link(x))
    cat("\n Mean formula:     ", formula_string(x$mformula))
    cat("\n Rvc formula:     ", formula_string(x$rformula))
    cat("\n observations:", nobs(x))
	cat("\n fixed effects: ", ncol(x$spec$mf$X))
    cat("\n------\n")

	mat <- cbind(coef(x),se(x))
	colnames(mat) <- c("Median","MAD")

	.printfr(mat,digits)

    mat <- as.matrix(x$rvcfit) # don't used as.matrix.rvcreg method b/c want access to mean_PPD

	smooth_nms <- grep("^smooth_precision\\[", colnames(mat), value = TRUE)
	smooth_mat <- mat[,smooth_nms,drop=F]
	smooth_ests <- .median_and_madsd(smooth_mat)

      cat("\nSmoothing terms:\n")
	.printfr(smooth_ests,digits)

}



# internal ----------------------------------------------------------------

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
