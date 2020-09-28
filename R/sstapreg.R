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

#' Create a sstapreg object
#'
#' The returned model object from the \pkg{rsstap} functions - methods can be 
#' called on this to identify STAP effects and any problems with model convergence 
#' or other diagnostics.
#' 
#' @param object A list provided by one of the \code{sstap_*} modeling functions.
#' @return A sstapreg object
#' @importFrom stats median mad
#'
sstapreg <- function(object){

	stapfit <- object$stapfit
	stanmat <- as.matrix(stapfit)
	stap_terms <- object$specification$term
	ynames <- grep("yhat",colnames(stanmat))
	nms <- colnames(object$mf$X)
	nms <- union(nms,Reduce(union,lapply(stap_terms,function(x) grep(paste0("^(s\\(|t2\\()",x),colnames(stanmat),value=T))))
	coefs <- apply(stanmat[,nms],2,median)
	covmat <- cov(stanmat[,nms])
	colnames(covmat) <-nms 
	ses <- apply(stanmat[,nms],2,mad)
	names(ses) <-nms 
	Nobs <- nrow((object$mf$X))
	
	y <- object$mf$y
	Z <- object$mf$X


    out <- list(
				coefficients = coefs, 
				ses = ses,
				fitted.values = apply(stanmat[,ynames],2,median),
				covmat = covmat,
				model = list(y=y,
							 Z=Z),
				specification = object$specification,
				formula = object$formula, 
				family = object$family,
				stapfit = stapfit,
				rstan_version = utils::packageVersion("rstan"),
				call = object$call, 
				stan_function = object$stan_function
		  )

	if(!is.null(object$glmod))
		out$glmod <- object$glmod

    structure(out, class = c("sstapreg"))

}
