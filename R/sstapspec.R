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
#' @param benvo Built Environment object - \code{\link[rbenvo]{benvo}} - containing data for model 
#' @return \code{\link{sstapspec}} object
#'
get_sstapspec <- function(f,benvo){

    with_bars <- lme4::findbars(f)
    f <- lme4::nobars(f)
	
	stapinfo <- get_bef_info(c("stap","sap","tap","gap"),f)
	if(is.null(stapinfo))
		stop("No covariates designated as ",paste0(c("stap","sap","tap","gap"),collapse=","),call. = F)

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
		   sstapspec(stapless_formula = as.formula(new_f, env = environment(f)),
					  fake_formula = fake_formula,
					  stap_mat = stapinfo,
					  benvo = benvo
			   )
		   )
}

get_bef_info <- function(nm_v, f){

	get_ics <- function(nm,f){
		which(all.names(f) %in% nm)
	}
	get_k <- function(nm,f){
		sig <- paste0("(^",nm,"\\()",collapse="|")
		tm <- attr(terms(f),"term.labels")
		strings <- tm[stringr::str_detect(tm,sig)]
		out <- stringr::str_extract(strings,", ?k ?= ?[1-9][0-9]? *\\)")
		out <- sapply(out,function(x) if(is.na(x)) return(")") else x)
		return(out)
	}
	get_names <- function(ics,f){
		all.names(f)[ics +1] 
	}
	get_indicator <- function(ics,nm,f){
		(all.names(f)[ics] %in% nm)*1
	}
	get_one_info <- function(nm,f){
	  nm <- c(nm,paste0(nm,"_bw"))
	  ics <- get_ics(nm,f)
	  if(length(ics)==0)
	    return(NULL)
	  nms <- get_names(ics,f)
	  bwi <- get_indicator(ics,nm[2],f)
	  k <- get_k(nm,f)
	  type <- nm[1]
	  out <- cbind(type,nms,k,bwi)
	  if(nrow(out)>0)
	    rownames(out) <-NULL
	  return(out)
	}
	out <- purrr::reduce(lapply(nm_v,function(x) get_one_info(x,f)),rbind)
	
	return(out)
}


#' Create STAP data structure
#' 
#' @export
#' @param stapless_formula from \code{\link{get_sstapspec}}
#' @param fake_formula list of ``fake'' formulas from \code{\link{get_sstapspec}}
#' @param stap_mat matrix of stap specification properties 
#' @param benvo Built Environment object - \code{\link[rbenvo]{benvo}} - containing data for model 
#'
sstapspec <- function(stapless_formula,fake_formula,stap_mat,benvo){

  component <- sapply(stap_mat[,1],function(x) switch(x,"sap"="Distance",
                                                          "tap"="Time",
                                                          "stap" = "Distance-Time",
                                                       "gap" = "Distance"))
	term <- stap_mat[,2]
	dimension <- as.integer(stringr::str_replace_na(stringr::str_extract(stap_mat[,3],"[1-9][0-9]?"),"-1"))
	between_within <- as.integer(stap_mat[,4])
	if(!(all(unique(term)==term)))
		stop("Only one BEF name may be assigned to a stap term e.g. no sap(foo) + tap(foo)\n
			 If you wish to model components this way create a different name e.g. sap(foo) + tap(foo_bar)")
	if(!all(term %in% rbenvo::bef_names(benvo)))
		stop("All stap terms must have data with corresponding name in benvo")


	jd <- purrr::pmap(list(term,component,fake_formula),
	                  function(x,y,z) {
						  temp_df <- rbenvo::joinvo(benvo,x,y,NA_to_zero = TRUE)
						  temp_df$tempix_ <- 1:nrow(temp_df)
	                    out <- mgcv::jagam(formula = z, 
	                                       family = gaussian(), 
	                                       data = temp_df,
										   file = tempfile(fileext = ".jags"), 
	                                       offset = NULL,
	                                       centred = FALSE,
	                                       diagonalize = FALSE)
	                    out$name <- x
						ix <- which(rbenvo::bef_names(benvo)==x)
						ranges <- list()
						if(y=="Distance"|y=="Distance-Time")
							ranges$Distance <- range(benvo$sub_bef_data[[ix]]$Distance,na.rm=T)
						if(y=="Time"|y=="Distance-Time")
							ranges$Time = range(benvo$sub_bef_data[[ix]]$Time,na.rm=T)
						out$ranges <- ranges
	                    return(out)
	                    })


	X <- purrr::pmap(list(term,component,between_within,jd),function(x,y,z,jdi){
			X <- create_X(x,y,z,
			              jdi$jags.data$X,benvo,
			              jdi$pregam$term.names)})

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
## Internal ------------


create_S_nms <- function(term,component){

	first <- sapply(component,function(x) 
	switch(x,"Distance"="s(",
		   "Time"="s(",
		   "Distance-Time"="t2("))
	out <- paste0(first,term,")")
	return(out)

}

get_new_f <- function(f,with_bars,not_needed){

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
	return(new_f)

}
