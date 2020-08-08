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
#' @param ... arguments for stan sampler
#' 
sstap_lm <- function(formula,
					  benvo,
					   ...){
  
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
	                          family = gaussian(),
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
