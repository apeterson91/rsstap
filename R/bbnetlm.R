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

#' Built Environment Networks Linear Model 
#'
#'
#' @export
#'
#' @param formula Similar as for \code{\link[stats]{lm}}. 
#' @param data  same as for \code{\link[stats]{lm}}
#' @param D sparse basis distance matrix
#' @param T translation matrix
#' @param iter number of iterations
#' @param warm_up number of warm up iterations
#' @param chains number of chains
#' @param seed random number generator
#'
bbnet_lm <- function(formula,
                   data,
                   D,
                   T_mat,
                   iter = 2E3,
                   warm_up = 1E3,
                   chains = 1L,
                   seed = NULL){

	if(is.null(seed))
		seed <- 1431
	
	call <- match.call(expand.dots = TRUE)
	mf <- match.call(expand.dots = FALSE)
	mf$formula <- formula
	m <- match(c("formula"),table = names(mf), nomatch=0L)
	mf <- mf[c(1L,m)]
	mf$data <- data
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf,parent.frame())
	mt <- attr(mf,"terms")
	if(is.empty.model(mt))
		stop("No intercept or predictors specified.",.call = FALSE)

	y <- model.response(mf,type = "any")
	Z <- model.matrix(mt,mf)
	if(attr(mt,"intercept")){
	  Z <- Z[,2:ncol(Z),drop=F]
	  has_intercept <- TRUE
	}
	else
	  has_intercept <- FALSE

	out <- bbnet_lm_fit(y = y,
	                   Z = Z, 
	                   DD = D,
	                   rpc = c(0,has_intercept,0,0,0,0),
	                   prior_means = c(0,0,0,0,0,0),
	                   prior_scales = rep(0,6),
	                   iter_max = iter,
	                   max_treedepth = 10, 
	                   warm_up = warm_up,
	                   seed = seed
	                   )

	
}
