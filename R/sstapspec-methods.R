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


#' \code{sstapspec} Methods
#'
#' @name sstapspec-methods
#' @aliases get_k get_terms_ix get_component has_bw get_smooth_obj 
#'
#' @param x sstapspec object
#' @template args-term
#' @details Methods for sstapspec objects. Function names are typically self-explanatory.


#' @rdname sstapspec-methods
#' @export
get_k <- function(x,term = NULL)
	UseMethod("get_k")


#'
#' 
#' @describeIn get_k retrieve dimension
#' @export
#'
get_k.sstapspec <- function(x, term = NULL){

	if(is.null(term))
		return(x$dimension)
	else
		return(x$dimension[get_terms_ix(x,term)])
}

#' @rdname sstapspec-methods
#' @export
get_terms_ix <- function(x,term= NULL)
	UseMethod("get_terms_ix")

#'
#' @describeIn get_terms_ix get term indices
#' @export
#'
get_terms_ix.sstapspec  <- function(x,term = NULL){
	if(is.null(term))
		return(1)
	else
		which(x$term==term)
}

#' @rdname sstapspec-methods
#' @export 
#'
get_component <- function(x,term = NULL)
	UseMethod('get_component')

#'
#' @describeIn get_component retrieve space-time designation 
#' @export
#'
get_component.sstapspec <- function(x, term = NULL){

	if(is.null(term))
		return(x$component)
	else
		return(x$component[get_terms_ix(x,term)])
}

#' @rdname sstapspec-methods
#' 
#' @export
#' 
has_bw <- function(x,term = NULL)
	UseMethod("has_bw")

#'
#' @describeIn has_bw has between-within decomposition
#' @export
#'
has_bw.sstapspec <- function(x,term = NULL){

	x$between_within[get_terms_ix(x,term)]
}


# Returns labels for given stap term, and vector of colnames
#
get_lbls <- function(x,term,colnms)
	UseMethod("get_lbls")

# 
get_lbls.sstapspec <- function(x,term,colnms){

	comp <- get_component(x,term)
	bw <- has_bw(x,term)
	patrn <- switch(comp,
		"Distance-Time"= "t2",
		"Distance" = "s",
		"Time" = "s"
	)
	if(bw)
		term <- paste0(term,"_(bw|wi)")
	patrn <- paste0(patrn,"\\(",term,"\\)\\.([1-9]$|[1-9][0-9]$)")
	out <- stringr::str_extract(colnms,patrn)
	out <- out[which(!is.na(out))]
	return(out)
}

has_any_staps <- function(x)
	UseMethod("has_any_staps")

has_any_staps <-function(x){
	return(any(x$component=="Distance-Time"))
}

get_range <- function(x,term,component)
	UseMethod("get_range")

get_range.sstapspec <- function(x,term,component){
	return(x$range[[which(x$term==term)]][[component]])
}

# Get Stap Grid
#
get_grid <- function(x,term, component = NULL)
	UseMethod("get_grid")

#
#
get_grid.sstapspec <- function(x,term,component = NULL){

	if(is.null(component))
		component <- get_component(x,term)

	if(component=="Distance-Time"){
		range_dist <- get_range(x,term,"Distance")
		range_time <- get_range(x,term,"Time")
	}else{
		range_ <- get_range(x,term,component)
	}

	generate_st_grid <- function(min_dist,max_dist,min_time,max_time){
		gd <- expand.grid(Distance = seq(from = min_dist, to = max_dist, by = 0.01),
						  Time = seq(from = min_time, to = max_time, by = 0.01)
		)
		return(as.data.frame(gd))
	}


	gd <- switch(component,
		   "Distance"= data.frame(Distance =seq(from = floor(range_[1]),
											   to = ceiling(range_[2]),
											   by = 0.01)),
		   "Distance-Time"= generate_st_grid(floor(range_dist[1]),ceiling(range_dist[2]),floor(range_time[1]),ceiling(range_time[2])),
		   "Time"= data.frame(Time = seq(from = floor(range_[1]),
											   to = ceiling(range_[2]),
											   by = 0.01))
	)

	return(gd)
}


#'
#' @rdname sstapspec-methods
#'
#' @export
#'
get_smooth_obj <- function(x,term = NULL)
	UseMethod("get_smooth_ob")

#' 
#' 
#' @describeIn get_smooth_ob Retrieve spline object
#' @export
#'
get_smooth_obj <- function(x,term = NULL){

	ix <- get_terms_ix(x,term)
	return(x$smooth_objs[[ix]])
}


get_stap <- function(x,term,component,beta,family)
	UseMethod("get_stap")

get_stap.sstapspec <- function(x,term,component,beta,family){

	gd <- get_grid(x,term,component)

	sob <- get_smooth_obj(x,term)

	mat <- mgcv::Predict.matrix(sob,gd)

	lbls <- get_lbls(x,term,colnames(beta))

	if(has_bw(x,term)){
		lbls_bw <- lbls[stringr::str_detect(lbls,"_bw")]
		lbls_wi <- lbls[stringr::str_detect(lbls,"_wi")]
		beta_ <- beta[,lbls_bw]
		eta_within = family$linkinv(tcrossprod(mat,beta_))
		beta_ <- beta[,lbls_wi] 
		eta <- family$linkinv(tcrossprod(mat,beta_))
		return(list(grid=gd,eta=eta,eta_within=eta_within))
	}

	beta <- beta[,lbls] 
	eta <- tcrossprod(mat,beta)

	return(list(grid=gd,
	            eta = eta))
}
