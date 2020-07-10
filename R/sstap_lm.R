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
#' @param weights vector of positive integer weights 
#' @param ... arguments for stan sampler
#' 
sstap_lm <- function(formula,
					 benvo,
					 weights = NULL,...){
  
	foo <- get_stapless_formula(formula)
	f <- foo$stapless_formula
	call <- match.call(expand.dots = TRUE)
	if(benvo@longitudinal){
		warning("This Benvo was constructed with longitudinal data \n but sstap_lm does not adjust for within subject correlation. \n Be advised that parameter standard errors may be overoptimistic.")
	}
	mf <- rbenvo::subject_design(benvo,f)
	stap_terms <- foo$stap_mat[,1]
	stap_components <- foo$stap_mat[,2]
	
	jd <- purrr::pmap(list(stap_terms,stap_components,foo$fake_formula),
	                  function(x,y,z) {
	                    out <- mgcv::jagam(formula = z, family = gaussian(), 
	                                       data = rbenvo::joinvo(benvo,x,y,
	                                                             NA_to_zero = TRUE), 
	                                       file = tempfile(fileext = ".jags"), 
	                                       weights = NULL, 
	                                       offset = NULL,
	                                       centred = FALSE,
	                                       diagonalize = FALSE)
	                    out$name <- x
						ix <- which(benvo@bef_names==x)
						ranges <- list()
						if(y=="Distance"|y=="Distance-Time")
							ranges$Distance <- range(benvo@bef_data[[ix]]$Distance,na.rm=T)
						if(y=="Time"|y=="Distance-Time")
							ranges$Time = range(benvo@bef_data[[ix]]$Time,na.rm=T)
						out$ranges <- ranges
	                    return(out)
	                    })
	
	X <- lapply(1:length(jd),function(i){
					X <- create_X(stap_terms[i],stap_components[i],jd[[i]]$jags.data$X,benvo,jd[[1]]$pregam$term.names)
					return(X)
			})
	## always only one smooth because jagam called for each stap term
	S <- lapply(1:length(jd),function(i){
			s <- jd[[i]]$jags.data$S1 
			if(!is.null(s))
				return(s)
			else
				return(diag(ncol(jd[[i]]$jags.data$X)))
			})
	

	
	

	sstapfit <- sstap_glm.fit(
	                          y = mf$y, 
	                          Z = mf$X,
	                          X = X,
	                          S = S,
	                          w = weights
	                          )
	
	fit <- sstapreg(
					list(stapfit = sstapfit$fit,
						 mf = mf,
						 smooths = lapply(jd,function(x) x$pregam$smooth),
						 ranges = lapply(jd,function(x) x$ranges),
						 weights = weights,
						 benvo = benvo,
						 Xs = lapply(jd,function(x) x$jags.data$X),
						 stap_terms = stap_terms,
						 stap_components = stap_components,
						 ind = sstapfit$ind,
						 call = call
					 )
					)

	return(fit)

}
