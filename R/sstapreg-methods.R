#' Methods for sstapreg objects
#' 
#' The methods documented on this page are actually some of the least important 
#' methods defined for \link[=sstapreg]{sstapreg} objects. The most 
#' important methods are documented separately, each with its own page. Links to
#' those pages are provided in the \strong{See Also} section, below.
#' 
#' @name sstapreg-methods
#' @aliases VarCorr ngrps sigma nsamples ranef
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
#' 
se <- function(object, ...) UseMethod("se")

#' @rdname sstapreg-methods
#' @export
se.sstapreg <- function(object, ...) {
  object$ses
}


#' @rdname sstapreg-methods
#'
#' @export
nobs.sstapreg <- function(object, ...){
	length(object$fitted.values)
}

#' @rdname sstapreg-methods
#'
#' @param x sstapreg object
#' @export
#' @importFrom nlme VarCorr
#'
VarCorr.sstapreg <- function(x){

	mat <- as.matrix(x)
	cnms <- .cnms(x)
	useSc <- "sigma" %in% colnames(mat)
    if (useSc) sc <- mat[,"sigma"] else sc <- 1
	Sigma <- colMeans(mat[,grepl("^Sigma\\[", colnames(mat)), drop = FALSE])
	nc <- vapply(cnms, FUN = length, FUN.VALUE = 1L)
	nms <- names(cnms)
	ncseq <- seq_along(nc)
	if (length(Sigma) == sum(nc * nc)) { # stanfit contains all Sigma entries
	spt <- split(Sigma, rep.int(ncseq, nc * nc))
	ans <- lapply(ncseq, function(i) {
	  Sigma <- matrix(0, nc[i], nc[i])
	  Sigma[,] <- spt[[i]]
	  rownames(Sigma) <- colnames(Sigma) <- cnms[[i]]
	  stddev <- sqrt(diag(Sigma))
	  corr <- cov2cor(Sigma)
	  structure(Sigma, stddev = stddev, correlation = corr)
	})       
	} else {  
	spt <- split(Sigma, rep.int(ncseq, (nc * (nc + 1)) / 2))
	ans <- lapply(ncseq, function(i) {
	  Sigma <- matrix(0, nc[i], nc[i])
	  Sigma[lower.tri(Sigma, diag = TRUE)] <- spt[[i]]
	  Sigma <- Sigma + t(Sigma)
	  diag(Sigma) <- diag(Sigma) / 2
	  rownames(Sigma) <- colnames(Sigma) <- cnms[[i]]
	  stddev <- sqrt(diag(Sigma))
	  corr <- cov2cor(Sigma)
	  structure(Sigma, stddev = stddev, correlation = corr)
	})    
	}
	names(ans) <- nms
	structure(ans, sc = mean(sc), useSc = useSc, class = "VarCorr.merMod")

}

#'
#' @rdname sstapreg-methods
#'
#' @export
#' @export ngrps
#' @importFrom lme4 ngrps
#' 
ngrps.sstapreg <- function(object, ...) {
  vapply(.flist(object), nlevels, 1)  
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
#' @export ranef
#' @importFrom lme4 ranef
#' @param benvo benvo object used to fit the model
#'
ranef.sstapreg <- function(object,benvo,...){
	.glmer_check(object)
	b_nms <- paste0("b[",make_b_nms(object$glmod$reTrms),"]")
	point_estimates <- summary(object)[b_nms,"50%"]
	  out <- ranef_template(object,benvo)
	  group_vars <- names(out)
  for (j in seq_along(out)) {
    tmp <- out[[j]]
    pars <- colnames(tmp) 
    levs <- rownames(tmp)
    levs <- gsub(" ", "_", levs) 
    for (p in seq_along(pars)) {
      stan_pars <- paste0("b[", pars[p], " ", group_vars[j],  ":", levs, "]")
      tmp[[pars[p]]] <- unname(point_estimates[stan_pars])
    }
    out[[j]] <- tmp
  }
  out
}

# Call lme4 to get the right structure for ranef objects
#' @importFrom lme4 lmerControl glmerControl lmer glmer 
ranef_template <- function(object,benvo) {
  
    new_formula <- object$specification$stapless_formula 
	if(family(object)$family =="gaussian")
		lme4_fun <- "lmer"
	else
		lme4_fun <- "glmer"

	cntrl_args <- list(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1))
	cntrl_args$check.conv.grad <- "ignore"
	cntrl_args$check.conv.singular <- "ignore"
	cntrl_args$check.conv.hess <- "ignore"
	cntrl_args$check.nlev.gtreq.5 <- "ignore"
	cntrl_args$check.nobs.vs.rankZ <- "ignore"
	cntrl_args$check.nobs.vs.nlev <- "ignore"
	cntrl_args$check.nobs.vs.nRE <- "ignore"
	if (lme4_fun == "glmer") {
      cntrl_args$check.response.not.const <- "ignore"
    }
  
  cntrl <- do.call(paste0(lme4_fun, "Control"), cntrl_args)
  
  fit_args <- list(
    formula = new_formula,
    data = benvo$subject_data,
    control = cntrl
  )
  
  if(family(object)$family != "gaussian"){
		fit_args$family <-family(object)
  }

	lme4_fit <- suppressWarnings(do.call(lme4_fun, args = fit_args))
	ranef(lme4_fit)
}

#'
#' @rdname sstapreg-methods
#' @export
#' @importFrom rstantools nsamples
nsamples.sstapreg <- function(object, ...) {
  posterior_sample_size(object)
}


#' Calculate WAIC
#' 
#' @export
#' @param x sstapreg object
#'
waic <- function(x)
	UseMethod("waic")

#' Calculate WAIC
#'
#' @describeIn waic calculate waic
#' @export
#' @importFrom stats dbinom dnorm dpois
#'
waic.sstapreg <- function(x){

	yhatmat <- as.matrix(x$stapfit)
	yhats <- grep("yhat",colnames(yhatmat))
	yhatmat <- t(yhatmat[,yhats])
	
	if(x$family$family=="binomial"){
	  if(is.matrix(x$model$y)){
		  nt <- rowSums(x$model$y)
		  y_ <- x$model$y[,1]
		  out <- LaplacesDemon::WAIC(t(apply(yhatmat,2,function(z) dbinom(x = y_, size = nt, prob = z,log = TRUE ))))
	  }else{
	    nt <- rep(1,length(x$model$y))
	    y_ <- x$model$y
	    out <- LaplacesDemon::WAIC(t(apply(yhatmat,2,function(z) dbinom(x = y_, size = nt, prob = z,log = TRUE ))))
	  }
	}else if(x$family$family=="gaussian"){
		sig <- as.matrix(x$stapfit)[,'sigma']
		out <- LaplacesDemon::WAIC(t(apply(yhatmat,2,function(z) dnorm(x = x$model$y,mean = z,sd = sig,log = TRUE ))))
	}else{
		##poisson
		out <- LaplacesDemon::WAIC(t(apply(yhatmat,2,function(z) dpois(x = x$model$y,lambda = z, log = TRUE))))
	}
	return(out)
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

.cnms <- function(object, ...) UseMethod(".cnms")
.cnms.sstapreg <- function(object, ...) {
  .glmer_check(object)
  object$glmod$reTrms$cnms
}

.flist <- function(object, ...) UseMethod(".flist")
.flist.sstapreg <- function(object, ...) {
  .glmer_check(object)
  as.list(object$glmod$reTrms$flist)
}

.glmer_check <- function(object) {
  if (!is.mer(object))
    stop("This method is for sstap_(g)lmer models only.", 
         call. = FALSE)
}
