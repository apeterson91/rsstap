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

#' Spline Spatial Temporal Aggregated Regression Linear Model Fit 
#'
#'
#' @export
#'
#' @param y vector of outcomes
#' @param Z matrix of subject level covariates
#' @param X list of Smooth design matrices
#' @param S list of Smooth penalty/precision matrices
#' @param ... arguments for stan sampler
#' 
sstap_glm.fit <- function(y,
						  Z,
						  X,
						  S,
						  w = NULL,
						  ...){
  
	N <- length(y)
	ncol_Z <- ncol(Z)
	
	K_smooth <- max(purrr::map_dbl(S,ncol))
	num_stap <- length(X)
	stap_penalties <- purrr::map_dbl(X,function(x) length(x))
	num_stap_penalties <- sum(stap_penalties)
	X <- purrr::map(X,function(x) do.call(cbind,x))
	stap_lengths <- purrr::map_dbl(X,function(x) ncol(x))
	ncol_smooth <- sum(stap_lengths)
	pen_ix <- matrix(0,nrow=num_stap_penalties,ncol=2)
	beta_ix <- matrix(0,nrow=num_stap,ncol=2)
	stap_pen_map <- matrix(0,nrow=num_stap,ncol=num_stap_penalties)
	startb <- 1
	startp <- 1
	cntr <- 1
	for(i in 1:num_stap){
	  end <- startb + stap_lengths[i] - 1L
	  beta_ix[i,1] <- startb
	  beta_ix[i,2] <- end
	  startb <- beta_ix[1,2] + 1L
	  for(j in 1:num_stap_penalties[i]){
	    stap_pen_map[i,j] <- cntr
	    cntr <- cntr + 1
	    ncol_S <- ncol(S[[i]])/num_stap_penalties[i]
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
	X <- cbind(Z,Reduce(cbind,X))
	qrc <- qr(X)
	Q <- qr.Q(qrc)
	R <- qr.R(qrc)
	R_inv <- solve(R)
	P <- ncol(Q)

		

  
  standata <- list(N = N,
                   ncol_Z = ncol_Z,
                   ncol_smooth = ncol_smooth,
                   num_stap = num_stap,
                   stap_lengths = array(stap_lengths),
                   stap_penalties = array(stap_penalties),
                   num_stap_penalties = num_stap_penalties,
                   pen_ix = pen_ix,
                   beta_ix = beta_ix,
                   P = P,
                   y = y, 
                   Q = Q,
                   R_inv = R_inv)

  pars <- c(
			"delta",
			"sstap_beta",
			"sigma",
			"tau",
			"yhat"
		  )

  stanfit <- stanmodels$sstap_continuous

  sampling_args <- set_sampling_args(
						object = stanfit,
						control = list(adapt_delta = 0.8,
						               max_treedepth = 10),
						pars = pars,
						data = standata,
						show_messages = FALSE,
						save_warmup = FALSE,
						...
						) 


  fit <- do.call(sampling,sampling_args)

  ind <- lapply(1:nrow(beta_ix),function(x) beta_ix[i,1]:beta_ix[i,2] )
  return(list (fit = fit, ind = ind))

}
