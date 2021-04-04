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

#' Spline Spatial Temporal Aggregated Generalized Mixed Effects Regression Model 
#'
#'
#' @export
#'
#' @param formula Similar as for \code{\link[lme4]{glmer}}. 
#' @param benvo built environment object from the rbenvo package containing the relevant data
#' @param family One of \code{\link[stats]{family}}  currently gaussian, binomial and poisson are implimented with identity, logistic and  log links currently.
#' @param QR boolean denoting whether or not to perform a QR decomposition on the design matrix, note that this is an experimental feature and bugs are still being worked out.
#' @param weights for unequal variances
#' @param ... optional arguments for stan sampler
#' 
sstap_glmer <- function(formula,
					   benvo,
					   family = gaussian(),
					   QR = TRUE,
					   weights = NULL,
					   ...){

	call <- match.call(expand.dots = TRUE)
	spec <- get_sstapspec(formula,benvo)
	f <- spec$stapless_formula
	mf <- rbenvo::longitudinal_design(benvo,f)
	
	sstapfit <- sstap_glm.fit(
	                          y = mf$y, 
	                          Z = mf$X,
	                          X = spec$X,
	                          S = spec$S,
	                          family = family,
	                          group = mf$glmod$reTrms,
							  QR = QR,
							  weights = weights,
							  ...
	                          )

	fit <- sstapreg(
					list(stapfit = sstapfit,
						 mf = mf,
						 benvo = benvo,
						 specification = spec,
						 call = call,
						 glmod = mf$glmod,
						 formula = formula,
						 family = family
					 )
					)


	return(fit)
}
