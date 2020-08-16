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

#' Spline Spatial Temporal Aggregated Predictor General Linear Model 
#'
#'
#' @export
#'
#' @param formula Similar as for \code{\link[stats]{glm}} with the addition of \code{stap},\code{sap} \code{tap} terms as needed
#' @param benvo built environment \code{\link[rbenvo]{Benvo}} object from the \code{rbenvo} package containing the relevant subject-BEF data
#' @param weights vector of positive integer weights 
#' @param family One of \code{\link[stats]{family}}  currently gaussian, binomial and poisson are implimented with identity, logistic and  log links currently.
#' @param QR boolean denoting whether or not to perform a QR decomposition on the design matrix.
#' @param ... arguments for stan sampler
#' 
sstap_glm <- function(formula,
					  benvo,
					  weights = NULL,
					  family = gaussian(),
					  QR = TRUE,
					   ...){
  
	validate_family(family)
	spec <- get_sstapspec(formula,benvo)
	f <- spec$stapless_formula
	call <- match.call(expand.dots = TRUE)
	mf <- rbenvo::subject_design(benvo,f)
	check_for_longitudinal_benvo(benvo)
	

	sstapfit <- sstap_glm.fit(
	                          y = mf$y, 
	                          Z = mf$X,
	                          X = spec$X,
	                          S = spec$S,
							  family = family,
							  QR = QR,
	                          ...
	                          )
	
	fit <- sstapreg(
					list(stapfit = sstapfit,
						 mf = mf,
						 weights = weights,
						 benvo = benvo,
						 specification = spec,
						 call = call,
						 formula = formula,
						 family = family
					 )
					)

	return(fit)

}
