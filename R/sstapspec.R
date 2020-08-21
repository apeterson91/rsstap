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


#' Retrieve STAP model specification from STAP model formula
#'
#' Get sstapspec object which details the various components of the sstap model specification
#'
#' @export
#' @param f formula from \code{\link{sstap_glm}}, \code{\link{sstap_glmer}}
#' @param benvo Built Environment object - \code{\link[rbenvo]{Benvo}} - containing data for model 
#' @return \code{\link{sstapspec}} object
#'
get_sstapspec <- function(f,benvo){

    with_bars <- lme4::findbars(f)
    f <- lme4::nobars(f)
	get_ics <- function(f,vec_var){
		which(all.names(f) %in% vec_var)
	}
	get_k <- function(strings){
	  out <- stringr::str_extract(strings,", ?k ?= ?[1-9][1-9]? *\\)")
	  out <- sapply(out,function(x) if(is.na(x)) return(")") else x)
	  return(out)
	}
	get_names <- function(f,ics){
		all.names(f)[ics +1] 
	}
	get_indicator <- function(f,ics,vec_var){
		(all.names(f)[ics] %in% vec_var)*1
	}
    stap_ics <- get_ics(f, c("stap","stap_bw"))
    sap_ics <- get_ics(f,c("sap","sap_bw"))
    tap_ics <- get_ics(f,c("tap","tap_bw"))
    if(!length(stap_ics) & !length(sap_ics) & !length(tap_ics))
        stop("No covariates designated as 'stap','sap',or 'tap'  in formula", .call = F)
	stap_nms <- get_names(f,stap_ics)
	stap_bw <- get_indicator(f,stap_ics, c("stap_bw"))
	sap_nms <- get_names(f,sap_ics)
	sap_bw <- get_indicator(f,sap_ics,c("sap_bw"))
    tap_nms <- get_names(f,tap_ics)
	tap_bw <- get_indicator(f,tap_ics,c("tap_bw")) 
	tms <- attr(terms(f),"term.labels")
	stap_tms <- tms[stringr::str_detect(tms,"(^stap\\()|(^stap_bw\\()")]
	tap_tms <- tms[stringr::str_detect(tms,"^tap\\(|(^tap_bw\\()")]
	sap_tms <- tms[stringr::str_detect(tms,"^sap\\(|(^sap_bw\\()")]
	stap_k <- get_k(stap_tms)
	tap_k <- get_k(tap_tms)
	sap_k <- get_k(sap_tms)
	
	if(length(stap_nms)>0){
		stap_nms <- cbind(stap_nms,"Distance-Time",stap_bw,stap_k)
	}
	if(length(sap_nms)>0)
		sap_nms <- cbind(sap_nms,"Distance",sap_bw,sap_k)
	if(length(tap_nms)>0)
		tap_nms <- cbind(tap_nms,"Time",tap_bw,tap_k)

	stap_mat <-rbind(stap_nms,sap_nms,tap_nms)

    not_needed <- c(stap_nms,sap_nms,tap_nms)
    formula_components <- all.vars(f)[!(all.vars(f) %in% not_needed)]
    if(!attr(terms(f),"intercept"))
        formula_components <- c(formula_components,"0")
    if(grepl("cbind",all.names(f))[2]){
        new_f1 <- paste0("cbind(",formula_components[1],", ",formula_components[2], ")", " ~ ")
        ix <- 3
    }
    else{
        new_f1 <- paste0(formula_components[1],' ~ ')
        ix <- 2
    }

    new_f2 <- paste(formula_components[ix:length(formula_components)],collapse = "+")
    new_f <- paste0(new_f1,new_f2)
	if(length(with_bars)){
		mer_f <- paste0(lapply(with_bars,function(x) paste0("(",deparse(x),")")),collapse = " + ")
		new_f <- paste0(new_f," + ",mer_f)
	}

	str <- purrr::map2(stap_mat[,2],stap_mat[,4],function(x,y) {
	  switch(x,
	         "Distance-Time"= paste0("t2(Distance,Time,bs='ps'",y ),
	         "Distance" = paste0("s(Distance,bs='ps'",y),
	         "Time"= paste0("s(Time,bs='ps'",y)
	         )
	  })

	fake_formula <- purrr::map(str,function(x) as.formula(paste0("ID~ -1 + ",paste0(x,collapse="+"))))

    return(
		   sstapspec(stapless_formula = as.formula(new_f, env = environment(f)),
					  fake_formula = fake_formula,
					  stap_mat = stap_mat,
					  benvo = benvo
			   )
		   )
}


#' Create STAP data structure
#' 
#' @export
#' @param stapless_formula from \code{\link{get_sstapspec}}
#' @param fake_formula list of ``fake'' formulas from \code{\link{get_sstapspec}}
#' @param stap_mat matrix of stap specification properties 
#' @param benvo Built Environment object - \code{\link[rbenvo]{Benvo}} - containing data for model 
#'
sstapspec <- function(stapless_formula,fake_formula,stap_mat,benvo){


	term <- stap_mat[,1]
	component <- stap_mat[,2]
	between_within <- as.integer(stap_mat[,3])
	dimension <- as.integer(stringr::str_replace(stap_mat[,4],"\\)","-1"))
	if(!(all(unique(term)==term)))
		stop("Only one BEF name may be assigned to a stap term e.g. no sap(foo) + tap(foo)\n
			 If you wish to model components this way create a different name e.g. sap(foo) + tap(foo_bar)")
	if(!all(term %in% benvo@bef_names))
		stop("All stap terms must have data with corresponding name in benvo")

	jd <- purrr::pmap(list(term,component,fake_formula),
	                  function(x,y,z) {
	                    out <- mgcv::jagam(formula = z, family = gaussian(), 
	                                       data = rbenvo::joinvo(benvo,x,y,
	                                                             NA_to_zero = TRUE), 
	                                       file = tempfile(fileext = ".jags"), 
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


	X <- purrr::pmap(list(term,component,between_within,jd),function(x,y,z,jdi){
			X <- create_X(x,y,z,
			              jdi$jags.data$X,benvo,
			              jdi$pregam$term.names)
						})

	S <- purrr::pmap(list(1:length(jd),between_within),function(ix,bw_){create_S(jd[[ix]],bw_)})
	names(S) <- create_S_nms(term,component)

	combine_list_entries <- function(l){
		if(any(sapply(l,is.list)))
			l <- Reduce(c,l)
		return(l)
	}

	X <- combine_list_entries(X)
	S <- combine_list_entries(S)
	



	out <- list(stapless_formula = stapless_formula,
				fake_formula = fake_formula,
				term = term,
				component = component,
				between_within = between_within,
				dimension = dimension,
				ranges = lapply(jd,function(x) x$ranges),
				X = X,
				S = S,
				smooth_objs = Reduce(c,lapply(jd,function(x) x$pregam$smooth))
				)

	structure(out,class=c("sstapspec"))
}


create_S_nms <- function(term,component){

	first <- sapply(component,function(x) 
	switch(x,"Distance"="s(",
		   "Time"="s(",
		   "Distance-Time"="t2("))
	out <- paste0(first,term,")")
	return(out)
}
