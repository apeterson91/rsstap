# This software is part of the rsstap package
# Copyright (C) 2020 Adam Peterson
#
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

#' Spline Spatial Temporal Aggregated Regression Linear Model 
#'
#' 
#' @export
#'
#' @param formula Similar as for \code{\link[stats]{lm}} with the addition of \code{stap},\code{sap} \code{tap} terms as needed
#' @param benvo built environment object from the rbenvo package containing the relevant data
#' @param QR boolean denoting whether or not to perform a QR decomposition on the design matrix.
#' @param weights for unequal variances
#' @param ... arguments for stan sampler
#' 
sstap_lm <- function(formula,
					  benvo,
					  QR = TRUE,
					  weights = NULL,
					   ...){
  
	spec <- get_sstapspec(formula,benvo)
	f <- spec$stapless_formula
	call <- match.call(expand.dots = TRUE)

	mf <- rbenvo::subject_design(benvo,f)
	check_for_longitudinal_benvo(benvo)
	if(!is.null(weights)){
		stopifnot(length(w)==length(mf$y))
		y <- sqrt(weights) * mf$y
		Z <- diag(sqrt(weights)) %*%  mf$X
		X <- diag(sqrt(weights)) %*% spec$X
	}else{
		Z <- mf$X
		y <- mf$y
		X <- spec$X
	}
	S <- spec$S
	


	sstapfit <- sstap_glm.fit(
	                          y = y, 
	                          Z = Z,
	                          X = X,
	                          S = S,
	                          family = gaussian(),
							  QR = QR,
	                          ...
	                          )
	
	fit <- sstapreg(
					list(stapfit = sstapfit,
						 mf = mf,
						 benvo = benvo,
						 specification = spec,
						 call = call,
						 formula = formula,
						 family = gaussian()
					 )
					)

	return(fit)

}
