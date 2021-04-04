# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3 # of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Spline Spatial Temporal Aggregated Regression Generalized Linear Model Fit 
#'
#'
#' @export
#'
#' @param y vector/matrix of outcomes
#' @param Z matrix of subject level covariates
#' @param X list of Smooth design matrices
#' @param S list of Smooth penalty/precision matrices
#' @param family One of \code{\link[stats]{family}}  currently gaussian, binomial and poisson are implimented with identity, logistic and  log links currently.
#' @param group A list, possibly of length zero (the default), but otherwise
#'   having the structure of that produced by \code{\link[lme4]{mkReTrms}} to
#'   indicate the group-specific part of the model. 
#' @param QR boolean denoting whether or not to perform a QR decomposition on the design matrix, note that this is an experimental feature and bugs are still being worked out.
#' @param weights for unequal variances
#' @param ... arguments for stan sampler
#' @importFrom Matrix t
#' 
sstap_glm.fit <- function(y,
						  Z,
						  X,
						  S,
						  family,
						  group = list(),
						  QR = F,
						  weights = NULL,
						  ...){
  
	if(is.matrix(y))
		N <- nrow(y)
	else
		N <- length(y)
	ncol_Z <- ncol(Z)
	
	if(!is.null(names(S)))
		S_nms <- names(S)
	else
		S_nms <-rep("",length(S))
	K_smooth <- max(purrr::map_dbl(S,ncol))
	num_stap <- length(X)
	stap_penalties <- purrr::map2_dbl(X,S,function(x,s) ncol(s)/ncol(x))
	num_stap_penalties <- sum(stap_penalties)
	stap_lengths <- purrr::map_dbl(X,function(x) ncol(x))
	X <- do.call(cbind,X)
	ncol_smooth <- sum(stap_lengths)
	pen_ix <- matrix(0,nrow=num_stap_penalties,ncol=2)
	beta_ix <- matrix(0,nrow=num_stap,ncol=2)
	stap_pen_map <- matrix(0,nrow=num_stap,ncol=max(stap_penalties))
	startb <- 1
	cntr <- 1
	for(i in 1:num_stap){
	  startp <- 1
	  end <- startb + stap_lengths[i] - 1L
	  beta_ix[i,1] <- startb
	  beta_ix[i,2] <- end
	  startb <- beta_ix[1,2] + 1L
	  for(j in 1: stap_penalties[i]){
	    stap_pen_map[i,j] <- cntr
	    cntr <- cntr + 1
	    ncol_S <- ncol(S[[i]])/stap_penalties[i]
	    end <- startp + ncol_S - 1L
	    pen_ix[stap_pen_map[i,j],1] <- startp
	    pen_ix[stap_pen_map[i,j],2] <- end
	    startp <- end + 1
	  }
	}
	S <- Reduce(rbind,lapply(S,function(x) {
	  num_zero <- K_smooth - ncol(x)
	  out <- cbind(x,matrix(0,ncol=num_zero,nrow=nrow(x)))
	  return(out)
	})
	)
	X <- cbind(Z,X)
	if(QR){
		qrc <- qr(X)
		Q <- qr.Q(qrc)
		R <- qr.R(qrc)
		R_inv <- qr.solve(qrc,Q)
	}else{
		Q <- X
		R_inv <- diag(ncol(X))
	}
	P <- ncol(Q)
  
  standata <- list(N = N,
                   ncol_Z = ncol_Z,
                   ncol_smooth = ncol_smooth,
                   num_stap = num_stap,
                   stap_lengths = array(stap_lengths),
                   stap_penalties = array(stap_penalties),
                   num_stap_penalties = num_stap_penalties,
                   stap_pen_map = stap_pen_map,
                   pen_ix = pen_ix,
                   beta_ix = beta_ix,
                   P = P,
                   y = y, 
                   Q = Q,
                   R_inv = R_inv)
  if(is.matrix(y)){
	  standata$y <- y[,1]
	  standata$num_trials <- rowSums(y) 
  }

  if(length(group)){
	## Code in this section has been largely adapted from
	## the rstanarm package. See github.com/stan-dev/rstanarm for more.
    check_reTrms(group)                                 
	Z <- t(group$Zt)                                     
	group <-
		pad_reTrms(Ztlist = group$Ztlist,
				   cnms = group$cnms,
				   flist = group$flist)
	Z <- group$Z
	p <- sapply(group$cnms,FUN = length)
	l <- sapply(attr(group$flist,"assign"), function(i) nlevels(group$flist[[i]]))
	t <- length(l)
	b_nms <- make_b_nms(group)
	g_nms <- unlist(lapply(1:t, FUN = function(i) {
				paste(group$cnms[[i]], names(group$cnms[i]), sep = "|")
	}))
	standata$t <- t
	standata$p <- as.array(p)
	standata$l <- as.array(l)
	standata$q <- ncol(Z)
	standata$len_theta_L <- sum(choose(p,2),p)
	parts <- rstan::extract_sparse_parts(Z)
	standata$num_non_zero <- length(parts$w)
	standata$w <- parts$w
	standata$v <- parts$v 
	standata$u <- parts$u
	standata$shape <- as.array(maybe_broadcast(1,t))
	standata$scale <- as.array(maybe_broadcast(1, t))
	standata$len_concentration <- sum(p[p > 1])
	standata$concentration <-as.array(maybe_broadcast(1,sum(p[p > 1])))
    standata$len_regularization <- sum(p > 1)
	standata$regularization <- as.array(maybe_broadcast(1,sum(p>1)))
	standata$special_case <- all(sapply(group$cnms, FUN = function(x) {length(x) == 1 && x == "(Intercept)"}))
  }else{
    standata$num_non_zero <- 0L
	standata$w <- double(0)
	standata$v <- integer(0)
	standata$u <- integer(0)
    standata$t <- 0L
    standata$p <- integer(0)
    standata$l <- integer(0)
    standata$q <- 0L
    standata$len_theta_L <- 0L
    standata$special_case <- 0L
    standata$shape <- standata$scale <- standata$concentration <-
      standata$regularization <- rep(0, 0)
    standata$len_concentration <- 0L
    standata$len_regularization <- 0L
  }
  if(is.null(weights)){
	  standata$has_weights <- 0L 
	  standata$weights <- numeric()
  }else{
	  stopifnot(length(weights)==length(y))
	  stopifnot(all(is.numeric(weights)))
	  standata$has_weights <- 1L
	  standata$weights <- weights 
  }


  pars <- c(
			"delta_coef",
			"sstap_beta",
			if(family$family=="gaussian") "sigma",
			"tau",
			if(length(group)) "b",
		    if(standata$len_theta_L) "theta_L"
			)

  stanfit <- pick_stanmodel(family$family)
  if(family$family=="binomial" && !is.matrix(y)){
	  standata$num_trials <- rep(1,N)
  }


  sampling_args <- set_sampling_args(
						object = stanfit,
						control = list(adapt_delta = 0.85,
						               max_treedepth = 10),
						pars = pars,
						data = standata,
						show_messages = FALSE,
						save_warmup = FALSE,
						...
						) 


	fit <- do.call(sampling,sampling_args)

	if(standata$len_theta_L){
		thetas_ref <- rstan::extract(fit,
									 pars = "theta_L", 
									 inc_warmup = FALSE,
									 permuted = FALSE)
		cnms <- group$cnms
		nc <- sapply(cnms, FUN = length)
		nms <- names(cnms)
		Sigma <- apply(thetas_ref, 1:2, FUN = function(theta) {
						   Sigma <- lme4::mkVarCorr(sc = 1, cnms, nc, theta, nms)
						   unlist(sapply(Sigma, simplify = FALSE,
										 FUN = function(x) x[lower.tri(x,TRUE)]))
							  })
		l <- length(dim(Sigma))
		end <- utils::tail(dim(Sigma), 1L)
		shift <- grep("^theta_L", names(fit@sim$samples[[1]]))[1] - 1L
		if(l==3) for (chain in 1:end) for (param in 1:nrow(Sigma)) {
			fit@sim$samples[[chain]][[shift + param]] <- Sigma[param, , chain]
		}
		else for (chain in 1:end) {
			fit@sim$samples[[chain]][[shift + 1]] <- Sigma[, chain]
		}
		Sigma_nms <- lapply(cnms, FUN = function(grp) {
								nm <- outer(grp, grp, FUN = paste, sep = ",")
								nm[lower.tri (nm , diag = TRUE)]
							  })
		for(j in seq_along(Sigma_nms)){
			Sigma_nms[[j]] <- paste0(nms[j], ":", Sigma_nms[[j]])
		}
		Sigma_nms <- unlist(Sigma_nms)
	}

	new_names <- c(colnames(X),
					if(family$family=="gaussian") "sigma",
					Reduce(c,purrr::map2(1:length(stap_penalties),stap_penalties,function(x,y) paste0("smooth_precision[",S_nms[x],1:y,"]"))),
                   if (length(group) && length(group$flist)) c(paste0("b[", b_nms, "]")),
                   if (standata$len_theta_L) paste0("Sigma[", Sigma_nms, "]"),
				   "log-posterior"
					)

    fit@sim$fnames_oi <- new_names


  return(fit)

}

make_b_nms <- function(group, m = NULL, stub = "Long") {
  group_nms <- names(group$cnms)
  b_nms <- character()
  m_stub <- NULL
  for (i in seq_along(group$cnms)) {
    nm <- group_nms[i]
    nms_i <- paste(group$cnms[[i]], nm)
    levels(group$flist[[nm]]) <- gsub(" ", "_", levels(group$flist[[nm]]))
    if (length(nms_i) == 1) {
      b_nms <- c(b_nms, paste0(m_stub, nms_i, ":", levels(group$flist[[nm]])))
    } else {
      b_nms <- c(b_nms, c(t(sapply(paste0(m_stub, nms_i), paste0, ":", 
                                   levels(group$flist[[nm]])))))
    }
  }
  return(b_nms)  
}

# Add extra level _NEW_ to each group
# 
# @param Ztlist ranef indicator matrices
# @param cnms group$cnms
# @param flist group$flist
pad_reTrms <- function(Ztlist, cnms, flist) {
  stopifnot(is.list(Ztlist))
  l <- sapply(attr(flist, "assign"), function(i) nlevels(flist[[i]]))
  p <- sapply(cnms, FUN = length)
  n <- ncol(Ztlist[[1]])
  for (i in attr(flist, "assign")) {
    levels(flist[[i]]) <- c(gsub(" ", "_", levels(flist[[i]])), 
                            paste0("_NEW_", names(flist)[i]))
  }
  for (i in 1:length(p)) {
    Ztlist[[i]] <- rbind(Ztlist[[i]], Matrix::Matrix(0, nrow = p[i], ncol = n, sparse = TRUE))
  }
  Z <- t(do.call(rbind, args = Ztlist))
  return(list(Z = Z,cnms = cnms, flist = flist))
}

check_reTrms <- function(reTrms) {
  stopifnot(is.list(reTrms))
  nms <- names(reTrms$cnms)
  dupes <- duplicated(nms)
  for (i in which(dupes)) {
    original <- reTrms$cnms[[nms[i]]]
    dupe <- reTrms$cnms[[i]]
    overlap <- dupe %in% original
    if (any(overlap))
      stop("Similar to rstanarm, rsstap does not permit formulas with duplicate group-specific terms.\n", 
           "In this case ", nms[i], " is used as a grouping factor multiple times and\n",
           dupe[overlap], " is included multiple times.\n", 
           "Consider using || or -1 in your formulas to prevent this from happening.")
  }
  return(invisible(NULL))
}
