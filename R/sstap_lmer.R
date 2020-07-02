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

#' Spline Spatial Temporal Aggregated Regression Model for Longitudinal Data 
#'
#'
#' @export
#'
#' @param formula Similar as for \code{\link[lme4]{glmer}}. 
#' @param stap_terms list of BEF terms
#' @param benvo built environment object from the rbenvo package containing the relevant data
#' @param weights vector of positive integer weights 
#' @param ... arguments for stan sampler
#' 
sstap_lmer <- function(formula,
					 stap_terms,
					 benvo,
					 weights = NULL,
					 ...){

  call <- match.call(expand.dots = TRUE)
  mf <- rbenvo::subject_design(benvo,formula)

  S_Xs <- lapply(stap_terms,function(x) smooth_matrix(benvo,x))
	
	sstapfit <- sstap_glm.fit(y = mf$y, 
							 Z = mf$X,
							 S_Xs = S_Xs,
							 w = weights)
	
	fit <- sstapreg(
	                list(stapfit = sstapfit,
	                     mf = mf,
						 weights = weights,
	                     S_Xs = S_Xs,
	                     benvo = benvo,
	                     stap_terms = stap_terms,
	                     call = call
	                     )
	                )
	return(fit)

}
