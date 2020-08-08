#This software is part of the rsstap package
#Copyright (C) 2020 Adam Peterson
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' Extract the posterior sample
#' 
#' Provides samples of the posterior parameters in a
#' matrix data structure
#' 
#' @method as.matrix sstapreg
#' @export
#' @param x sstaprego object
#' @param pars vector of parameter names to subset from samples matrix
#' @param ... Ignored.
#'   
#' @return A matrix  the dimensions of which depend on
#'   \code{pars} and \code{regex_pars}..
#' 
#' @seealso \code{\link{sstapreg-methods}}
#' 
#' @examples
#' \donttest{
#' if (!exists("example_model")) example(example_model)
#' # Extract posterior sample after MCMC
#' draws <- as.matrix(example_model)
#' print(dim(draws))
#' 
#' # For example, we can see that the median of the draws for the intercept 
#' # is the same as the point estimate rstanarm uses
#' print(median(draws[, "(Intercept)"]))
#' print(example_model$coefficients[["(Intercept)"]])
#' }
#' 
as.matrix.sstapreg <- function(x, ..., pars = NULL) {
  
	mat <- as.matrix(x$stapfit)
	nms <- grep("yhat",colnames(mat),invert=T,value=T)
	if(!is.null(pars))
	  nms <- union(nms,pars)

	mat <- mat[, nms, drop = FALSE]
	return(mat)
}
