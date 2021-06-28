#  This software is part of the rsstap package
#  Copyright (C) 2020 Adam Peterson
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' Retrieve STAP-RVC model specification from STAP model formula
#'
#' Get rvcspec object which details the various components of the rvc model specification
#'
#' @export
#' @param mf formula for observation level mean estimation
#' @param rf formula for spatio-temporal level coefficient estimation
#' @param benvo Built Environment object - \code{\link[rbenvo]{benvo}} - containing data for model 
#' @return \code{\link{rvcspec}} object
#'
get_rvcspec <- function(mf,rf,benvo){

    with_bars <- lme4::findbars(mf)
    f <- lme4::nobars(mf)
	
	stapinfo <- get_bef_info(c("stap","sap","tap"),f)
	if(is.null(stapinfo))
		stop("No covariates designated as ",paste0(c("stap","sap","tap"),collapse=","),call. = F)

	new_f <- get_new_f(f,with_bars,stapinfo[,2])

	str <- purrr::map2(stapinfo[,1],stapinfo[,3],function(x,y) {
	  switch(x,
	         "stap"= paste0("t2(Distance,Time,bs='ps'",y ),
	         "sap" = paste0("s(Distance,bs='ps'",y),
	         "tap"= paste0("s(Time,bs='ps'",y),
			 "gap" = paste0("t2(Distance,Rank,bs='ps'",y)
	         )
	  }
	)

	fake_formula <- purrr::map(str,function(x) as.formula(paste0("tempix_ ~ -1 + ",paste0(x,collapse="+"))))

    return(
		   rvcspec(stapless_formula = as.formula(new_f, env = environment(f)),
				   fake_formula = fake_formula,
				   rf = rf,
				   components = get_components(stapinfo[,1]),
				   terms = stapinfo[,2],
				   dimensions = get_dimension(stapinfo[,3]),
				   between_withins <- as.integer(stapinfo[,4]),
				  benvo = benvo
			   )
		   )
}

get_components <- function(x){

	out <- sapply(x,
				function(y) switch(y,"sap"="Distance",
								   "tap"="Time",
								   "stap" = "Distance-Time",
								   "gap" = "Distance"))
	names(out) <- NULL
	return(out)
}


get_dimension <- function(x){

	out <- as.integer(stringr::str_replace_na(stringr::str_extract(x,"[1-9][0-9]?"),"-1"))
	return(out)
}



#' Random Varying Coefficients Model Specification
#'
#' 
#' @export
#' @param stapless_formula from \code{\link{get_rvcspec}}
#' @param fake_formula list of ``fake'' formulas from \code{\link{get_sstapspec}}
#' @param rf formula for spatio-temporal level coefficient estimation
#' @param components Distance/Time dimension to be modeled
#' @param terms terms of BEFs
#' @param dimensions dimension of spline basis expansion
#' @param between_within  indicator variable determining whether to perform decomposition
#' @param benvo built environment object - \code{\link[rbenvo]{benvo}} - containing data for model
#' 
rvcspec <- function(stapless_formula,fake_formula,rf,components,terms,dimensions,between_withins,benvo){


	jdf <- rbenvo::joinvo(benvo,terms,components, NA_to_zero=TRUE)
	ids <- rbenvo::get_id(benvo)
	D <- model.frame(rf,jdf %>% mutate_all(function(x) replace_na(x,0)))
	D <- model.matrix(D,jdf %>% mutate_all(function(x) replace_na(x,0)))
	if(any(colnames(D)=="(Intercept)"))
	  D <- D[,stringr::str_detect(colnames(D),"(Intercept)",negate=T),drop=F]


	jd <- mgcv::jagam(formula = fake_formula[[1]],
					  family = gaussian(),
					  data = jdf %>% dplyr::mutate(tempix_ = rnorm(dplyr::n())), 
					  file = tempfile(fileext = '.jags'),
					  offset = NULL,
					  centred = FALSE,
					  diagonalize = FALSE)

	jd$ranges <- range(jdf[,components,drop=T],na.rm=T)

	AggMat <- aggregation_matrix(benvo,terms,components)

	mf <- rbenvo::subject_design(benvo,stapless_formula)

	X <- jd$jags.data$X

	colnames(X) <- paste0("s(",terms,".",1:ncol(X),")")


	out <- list(
				stapless_formula = stapless_formula,
				fake_formula = fake_formula,
				term = terms,
				component = components,
				between_within = between_withins,
				dimension = dimensions,
				ranges = jd$ranges,
				AggMat = AggMat,
				X =  X,
				S1 = jd$jags.data$S1[1:10,1:10],
				S2 = jd$jags.data$S1[1:10,11:20],
				D = D,
				smooth_objs = jd$pregam$smooth,
				mf = mf
				)

	return(structure(out,class = c("rvcspec")))

}

#' Retrieves names for mean and rvc vector coefficients
#' 
#' 
get_coefnames <- function(x,type='all') UseMethod("get_coefnames")

get_coefnames.rvcspec <- function(x,type='all'){

	fixef <- colnames(x$mf$X)
	smooth <- colnames(x$X)
	rvc <- colnames(x$D)
	if(type=='all')
		return(c(fixef,smooth,rvc))
	if(type=='fixef')
		return(fixef)
	if(type=='smooth')
		return(smooth)
	if(type=="rvc")
		return(rvc)
}
