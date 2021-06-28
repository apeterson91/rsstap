#  This software is part of the rsstap package
#  Copyright (C) 2020 Adam Peterson
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Set arguments for sampling 
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
  args$control <- list(adapt_delta = 0.85,
					   max_treedepth = 10)

  return(args)
}



create_S <- function(jg,bw){

	if(is.null(jg$jags.data$S1)){
		## t2(Distance,Time)
		S <- Reduce(cbind,lapply(jg$pregam$S,function(x){
				   ix_nonzero <- which(diag(x)!=0)
				   out <- matrix(0,nrow = ncol(jg$jags.data$X),
								 ncol = ncol(jg$jags.data$X))
				   diag(out)[ix_nonzero] <- 1
				   return(out)
		   }))
	}else{
		##s(Distance) or s(Time)
		S <- jg$jags.data$S1
	}
	if(bw)
		return(list(S,S))

	return(S)
}


create_X <- function(stap_term,stap_component,calc_bw,raw_X,benvo,lbls = NULL){


	X <- rbenvo::aggrenvo(benvo,raw_X,stap_term,stap_component)
	if(!is.null(lbls)){
		if(stap_component == "Distance-Time")
		  stap_component_ <- "Distance,Time"
		else
		  stap_component_ <- stap_component
		lbls_ <- stringr::str_replace(lbls,stap_component_,stap_term)
		if(calc_bw)
			lbls_ <- lapply(c("_bw","_wi"), function(x) stringr::str_replace(lbls_,stap_term,paste0(stap_term,x)))
		else{
			colnames(X) <- lbls_
			return(X)
		}
	}
	if(calc_bw){
		X <- rbenvo::bwinvo(benvo,X)
		if(!is.null(lbls))
			X <- purrr::map2(lbls_,X,function(x,y){
					colnames(y) <- x
					return(y) 
			   })
	}

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
	if(rbenvo::is.longitudinal(benvo)){
		warning("This benvo was constructed with longitudinal data \n but sstap_lm does not adjust for within subject correlation. \n Be advised that parameter standard errors may be overoptimistic.")
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



# Maybe broadcast 
#
# @param x A vector or scalar.
# @param n Number of replications to possibly make. 
# @return If \code{x} has no length the \code{0} replicated \code{n} times is
#   returned. If \code{x} has length 1, the \code{x} replicated \code{n} times
#   is returned. Otherwise \code{x} itself is returned.
maybe_broadcast <- function(x, n) {
  if (!length(x)) {
    rep(0, times = n)
  } else if (length(x) == 1L) {
    rep(x, times = n)
  } else {
    x
  }
}

is.mer <- function(x){
	if(!is.null(x$glmod))
		return(TRUE)
	else
		return(FALSE)
}

create_unique_ID_mat <- function(id_one,id_two = NULL){
	tmp <- paste0(id_one,"_",id_two)
	lvls <- unique(tmp)
	new_id <- factor(tmp,levels=lvls)
	Matrix::fac2sparse(new_id)
}
