# This software is part of the rsstap package
# Copyright (C) 2020 Adam Peterson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Create a rvcreg object
#'
#' The returned model object from the \pkg{rsstap} functions - methods can be 
#' called on this to identify STAP effects and any problems with model convergence 
#' or other diagnostics.
#' 
#' @param object A list provided by one of the \code{sstap_*} modeling functions.
#' @return A rvcreg object
#' @importFrom stats median mad
#'
rvcreg <- function(object){

	stanmat <- as.matrix(object$rvcfit)

	nms <- get_coefnames(object$spec)
	coefs <- apply(stanmat[,nms],2,median)
	ses <- apply(stanmat[,nms],2,median)
	
    out <- list(
				stanmat = stanmat,
				coefficients = coefs,
				ses = ses,
			#	fitted.values = fitted.values,
				#covmat = covmat,
			#	model = list(y=y,
			#				 Z=Z),
				spec = object$spec,
				mformula = object$mean_formula, 
				rformula = object$rvc_formula,
				family = stats::gaussian(),
#				family = object$family,
				rvcfit = object$rvcfit,
				rstan_version = utils::packageVersion("rstan")
#				call = object$call, 
#				stan_function = object$stan_function
		  )

    structure(out, class = c("rvcreg"))

}
