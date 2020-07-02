#' Create a sstapreg object
#'
#' @param object A list provided by one of the \code{sstap_*} modeling functions.
#' @return A sstapreg object
#' @importFrom stats median mad
#'
sstapreg <- function(object){

	stapfit <- object$stapfit
	stanmat <- as.matrix(stapfit)
	del_names <- grep("delta",colnames(stanmat))
	ynames <- grep("yhat",colnames(stanmat))
	Z_coefs <- apply(stanmat[,del_names],2,median)
	Z_names <- colnames(object$mf$X)
	names(Z_coefs) <- Z_names
	covmat <- cov(stanmat[,del_names])
	colnames(covmat) <- Z_names
	stap_terms <- object$stap_terms
	ses <- apply(stanmat[,del_names],2,mad)
	names(ses) <- Z_names
	Nobs <- nrow((object$mf$X))
	
	if(!is.null(object$weights)){
		L <- sqrt(diag(object$weights))
		Z <- L %*% (object$mf$X)
		y <- L %*% object$mf$y
	}
	else{
		y <- object$mf$y
		Z <- object$mf$X
	}


    

    out <- list(
        coefficients = Z_coefs, 
        ses = ses,
        fitted.values = apply(stanmat[,ynames],1,median),
        covmat = covmat,
        model = list(y=y,
                     Z=Z,
					 smooths = object$smooths,
					 benvo = object$benvo),
        formula = object$formula, 
        stapfit = stapfit,
        stap_terms = stap_terms,
		ind = object$ind,
		ranges = object$ranges,
		stap_components = object$stap_components,
        rstan_version = utils::packageVersion("rstan"),
        call = object$call, 
        # sometimes 'call' is no good (e.g. if using do.call(stap_glm, args)) so
        # also include the name of the modeling function (for use when printing,
        # etc.)
        stan_function = object$stan_function
      )
	if(!is.null(object$weights))
		out$model$weights <- object$weights

    structure(out, class = c("sstapreg"))

}
