
# Set arguments for sampling 
#
# Prepare a list of arguments to use with \code{rstan::sampling} via
# \code{do.call}.
#
# @param object The sstapfit object to use for sampling.
# @param user_dots The contents of \code{...} from the user's call to
#   the \code{stan_*} modeling function.
# @param user_adapt_delta The value for \code{adapt_delta} specified by the
#   user.
# @param prior Prior distribution list (can be NULL).
# @param ... Other arguments to \code{\link[rstan]{sampling}} not coming from
#   \code{user_dots} (e.g. \code{data}, \code{pars}, \code{init}, etc.)
# @return A list of arguments to use for the \code{args} argument for 
#   \code{do.call(sampling, args)}.
set_sampling_args <- function(object, user_dots = list(), 
                               ...) {
  args <- list(object = object, ...)
  unms <- names(user_dots)
  for (j in seq_along(user_dots)) {
    args[[unms[j]]] <- user_dots[[j]]
  }
  args$control <- list(adapt_delta = 0.8,
					   max_treedepth = 10)

  return(args)
}

#' get_stapless_formula
#'
#' Get formula for typical covariates
#'
#' @param f formula from stap_glm
#' @return formula without  stap(),sap(),tap() components
#'
get_stapless_formula <- function(f){
    
    with_bars <- f
    f <- lme4::nobars(f)
	get_ics <- function(f,vec_var){
		which(all.names(f) %in% vec_var)
	}
    stap_ics <- get_ics(f, c("stap","stap_bw"))
    sap_ics <- get_ics(f,c("sap","sap_bw"))
    tap_ics <- get_ics(f,c("tap","tap_bw"))
    if(!length(stap_ics) & !length(sap_ics) & !length(tap_ics))
        stop("No covariates designated as 'stap','sap',or 'tap'  in formula", .call = F)
    stap_nms <- all.names(f)[stap_ics + 1]
	stap_bw <- (all.names(f)[stap_ics] %in% c("stap_bw"))*1
    sap_nms <- all.names(f)[sap_ics + 1]
	sap_bw <- (all.names(f)[sap_ics] %in% c("sap_bw"))*1
    tap_nms <- all.names(f)[tap_ics + 1]
	tap_bw <- (all.names(f)[tap_ics] %in% c("tap_bw"))*1
	if(length(stap_nms)>0){
		stap_nms <- cbind(stap_nms,"Distance-Time",stap_bw)
	}
	if(length(sap_nms)>0)
		sap_nms <- cbind(sap_nms,"Distance",sap_bw)
	if(length(tap_nms)>0)
		tap_nms <- cbind(tap_nms,"Time",tap_bw)

	stap_mat <-rbind(stap_nms,sap_nms,tap_nms)

    not_needed <- c(stap_nms,sap_nms,tap_nms)
    formula_components <- all.vars(f)[!(all.vars(f) %in% not_needed)]
    if(!attr(terms(f),"intercept"))
        formula_components <- c(formula_components,"0")
    if(grepl("cbind",all.names(f))[2]){
        new_f1 <- paste0("cbind(",formula_components[1],", ",formula_components[2], ")", " ~ ")
        ix <- 3
    }
    else{
        new_f1 <- paste0(formula_components[1],' ~ ')
        ix <- 2
    }

    new_f2 <- paste(formula_components[ix:length(formula_components)],collapse = "+")
    new_f <- paste0(new_f1,new_f2)

	str <- purrr::map(stap_mat[,2],function(x) {
	  switch(x,
	         "Distance-Time"= "t2(Distance,Time,bs='ps')", 
	         "Distance" = "s(Distance,bs='ps')",
	         "Time"= "s(Time,bs='ps')")
	  })

	fake_formula <- purrr::map(str,function(x) as.formula(paste0("ID~ -1 + ",paste0(x,collapse="+"))))

    return(
		   list(stapless_formula = as.formula(new_f, env = environment(f)),
				fake_formula = fake_formula,
				stap_mat = stap_mat
			   )
		   )
}


create_X <- function(stap_term,stap_component,calc_bw,raw_X,benvo,lbls = NULL){


	X <- rbenvo::aggrenvo(benvo,raw_X,stap_term,stap_component)
	if(!is.null(lbls))
		colnames(X) <- lbls 

	return(X)
}
# Get the posterior sample size
#
# @param x A stanreg object
# @return the posterior sample size (or size of sample from approximate posterior)
posterior_sample_size <- function(x) {
  pss <- x$stapfit@sim$n_save
  sum(pss - x$stapfit@sim$warmup2)
}


# Checks to see if benvo is longitudinal for _lm functions
# Prints warning message accordingly
#
# @param x benvo
# @return warning message if appropriate
check_for_longitudinal_benvo <- function(benvo){
	if(benvo@longitudinal){
		warning("This Benvo was constructed with longitudinal data \n but sstap_lm does not adjust for within subject correlation. \n Be advised that parameter standard errors may be overoptimistic.")
	}
}

# validates family is in list of reccomended families
validate_family <- function(family){
	if(!(family$family %in% c("gaussian","poisson","binomial")))
		stop("family must be in one of c('gaussian','poisson','binomial')")
	if(family$family =="gaussian" & family$link != "identity")
		warning("Currently only identity link supported for gaussian family")
	if(family$family =="binomial" & family$link != "logit")
		warning("Currently only logit link supported for binomial family")
	if(family$family =="poisson" & family$link != "log")
		warning("Currently only log link supported for binomial family")
}

pick_stanmodel <- function(family){
	if(family == "gaussian")
		return(stanmodels$sstap_continuous)
	else if(family == "binomial")
		return(stanmodels$sstap_binomial)
	else
		return(stanmodels$sstap_count)
}


