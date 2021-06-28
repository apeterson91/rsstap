

#' Plot Spline Spatial Temporal Aggregated Predictors
#' 
#' Plots the smooth curve on a grid over the range of the stap smooth predictor
#' If no stap term is selected, the first stap component is plotted by default.
#' @export
#' @param x sstapreg object
#' @param stap_term optional string for name of BEF smooth function to plot. 
#' Alternatively plots first BEF smooth function 
#' @param p probability mass contained within uncertainty interval
#' @param grid by default is NULL corresponding to a .1 step grid from the min-max of the covariate space
#' @param ... ignored
#' 
plot.sstapreg <- function(x,stap_term = NULL, p = 0.95, grid = NULL,...){

	# to pass R CMD Check
	Distance <- Time <- Median <- Parameters <- 
		Grid <- Lower <- Upper <-  . <-  .data <- NULL
	spec <- x$specification

	if(is.null(stap_term)){
		ix <- 1
		stap_term <- spec$term[ix]
		component <- get_component(spec,stap_term)
	}
	else if(stap_term %in% spec$term){
		ix <- which(spec$term==stap_term)
		component <- get_component(spec,stap_term)
	}
	else
		stop("stap_term must be NULL or one of the stap terms in the model object")

	beta <- as.matrix(x$stapfit)
	gd_eta <- get_stap(spec,stap_term,component,beta,x$family,grid)
	gd <- gd_eta$grid
	eta <- gd_eta$eta

	pltdf <- get_pltdf(gd,eta,p,component)

	if(component %in% c("Distance","Time")){
	  
		pltdf %>% ggplot2::ggplot(ggplot2::aes(x = .data[[component]], y = Median )) + 
			ggplot2::geom_line() + 
			ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),alpha=0.3) + 
			ggplot2::theme_bw() + 
			ggplot2::geom_hline(ggplot2::aes(yintercept = 0),
								linetype=2,color='red') + 
			ggplot2::xlab(component) + 
			ggplot2::ggtitle(stap_term) + 
			ggplot2::ylab("") -> pl
	}else{
	  pltdf %>% ggplot2::ggplot(ggplot2::aes(x=Distance,y=Time,z=Median)) + 
		  ggplot2::theme_bw() + 
	    ggplot2::geom_contour()  ->pl
	}

	if(has_bw(spec,stap_term)){

		pltdf2 <- get_pltdf(gd,gd_eta$eta_within,p,component) %>% 
			dplyr::mutate(Parameters = "Within")
		pltdf %>% dplyr::mutate(Parameters=="Between") %>% 
			rbind(.,pltdf2) -> pltdf

		if(component %in% c("Distance","Time")){
			pltdf %>% ggplot2::ggplot(ggplot2::aes(x=.data[[component]],y=Median)) + 
					  ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),alpha=0.3) + 
					  ggplot2::theme_bw() + 
					  ggplot2::theme(strip.background=ggplot2::element_blank()) + 
					  ggplot2::facet_wrap(~Parameters) -> pl2
		}else{
			pltdf %>% ggplot2::ggplot(ggplot2::aes(x=Distance,y=Time,z=Median)) + 
					  ggplot2::geom_contour() + 
					  ggplot2::theme_bw() + 
					  ggplot2::theme(strip.background=ggplot2::element_blank()) + 
					  ggplot2::facet_wrap(~Parameters) -> pl2
		}
	

	return(pl2)
	}
	
	return(pl)
}

#' Spatial-Temporal Effects DataFrame
#' 
#' Returns a dataframe with the betas evaluated at either a default or given grid points
#' If no stap term is selected, the first stap component is plotted by default.
#' @export
#' @param x sstapreg object
#' @param stap_term optional string for name of BEF smooth function to plot. 
#' Alternatively plots first BEF smooth function 
#' @param p probability mass contained within uncertainty interval
#' @param grid by default is NULL corresponding to a .1 step grid from the min-max of the covariate space
plotdf <- function(x,stap_term = NULL, p = 0.95, grid = NULL)
	UseMethod("plotdf")

#' @describeIn  plotdf
#' @export
plotdf.sstapreg <- function(x,stap_term = NULL, p = 0.95, grid = NULL){

	# to pass R CMD Check
	Distance <- Time <- Median <- Parameters <- 
		Grid <- Lower <- Upper <-  . <-  .data <- NULL
	spec <- x$specification

	if(is.null(stap_term)){
		ix <- 1
		stap_term <- spec$term[ix]
		component <- get_component(spec,stap_term)
	}
	else if(stap_term %in% spec$term){
		ix <- which(spec$term==stap_term)
		component <- get_component(spec,stap_term)
	}
	else
		stop("stap_term must be NULL or one of the stap terms in the model object")

	beta <- as.matrix(x$stapfit)
	gd_eta <- get_stap(spec,stap_term,component,beta,x$family,grid)
	gd <- gd_eta$grid
	eta <- gd_eta$eta

	pltdf <- get_pltdf(gd,eta,p,component,eta_wi = gd_eta$eta_wi)
	return(pltdf)
}




#' 3D plots for rsstap models
#'
#' @export
#' @param x sstapreg object
#' @param stap_term string argument for which \code{stap()} bef term to plot
#' @param grid  optional grid argument by default is NULL corresponding to a .1 step grid from the min-max of the covariate space
#'
plot3D <- function(x,stap_term = NULL,grid = NULL)
	UseMethod("plot3D")

#'
#' @describeIn plot3D 3d stap plot
#' @export 
#'
plot3D.sstapreg <- function(x,stap_term = NULL, grid = NULL){

	spec <- x$specification
	if(!has_any_staps(spec))
		stop("Model has no stap terms specified")

	if(is.null(stap_term)){
		ix <- which(spec$component=="Distance-Time")[1]
		stap_term <- spec$term[ix]
	}
	else if(stap_term %in% spec$term){
		## check term is a stap_term
		stopifnot(get_component(spec,stap_term)=="Distance-Time")
	}else{
		stop("Term specified is not included in model")
	}

	beta <- as.matrix(x$stapfit)

	gd_eta <- get_stap(spec,stap_term,"Distance-Time",beta,x$family,grid)
	gd <- gd_eta$grid
	eta <- gd_eta$eta

  	dplyr::tibble(Distance = gd$Distance,
				  Time = gd$Time,
  				  Lower = apply(eta,1,function(x) quantile(x,0.025)),
  				  Exposure = apply(eta,1,median),
  				  Upper = apply(eta,1,function(x) quantile(x,0.975))
				  ) -> pltdf
  
	pl <- plotly::plot_ly(pltdf,x = ~Distance, y = ~Time, z = ~Exposure,
						  mode = "markers", 
						  type ='scatter3d') 
	if(has_bw(spec,stap_term)){
		dplyr::tibble(Distance = gd$Distance,
					  Time = gd$Time,
					  Lower = apply(gd_eta$eta_within,1,function(x) quantile(x,0.025)),
					  Exposure = apply(gd_eta$eta_within,1,median),
					  Upper = apply(gd_eta$eta_within,1,function(x) quantile(x,0.975))
					  ) -> pltdf
	  
		pl2 <- plotly::plot_ly(pltdf,x = ~Distance, y = ~Time, z = ~Exposure,
							  mode = "markers", 
							  type ='scatter3d') 

		return((list(Between=pl,Within=pl2)))
	}

	return(pl)
}

#' Posterior Predictive Checks 
#'
#' @export
#' @keywords internal
#' @param x a sstapreg object
#' @param num_reps number of yhat samples to plot
#'
ppc <- function(x,num_reps = 20)
	UseMethod("ppc")

#' Posterior Predictive Checks
#'
#' @export
#' @describeIn ppc
#'
ppc.sstapreg <- function(x,num_reps = 20){

	Samples <- Parameter <- iteration_ix <- NULL

	samp <- sample(1:nsamples(x),size=num_reps)

	if(is.matrix(x$model$y))
	  y <- x$model$y[,1] / rowSums(x$model$y)
	else
	  y <- x$model$y

	yhatmat <- posterior_predict(x,'response')
	yhatmat <- yhatmat[samp,]
	colnames(yhatmat) <- paste0("yhat_",colnames(yhatmat))

	pltdf <- dplyr::as_tibble(yhatmat,silent=TRUE) %>% 
		dplyr::mutate(iteration_ix = 1:dplyr::n()) %>%
		tidyr::pivot_longer(dplyr::contains('yhat'),
							names_to = "Parameter",
							values_to="Samples") %>% 
						  dplyr::mutate(Parameter = 'yrep')
	
	pltdf <- rbind(pltdf,
				   dplyr::tibble(iteration_ix = 0, Parameter='y',Samples= y ))
	
	p <- pltdf %>% 
		ggplot2::ggplot(ggplot2::aes(x=Samples,color=Parameter,group=iteration_ix)) + 
		ggplot2::geom_density() + ggplot2::theme_bw()+ ggplot2::theme(legend.title=ggplot2::element_blank()) +   
		ggplot2::scale_colour_manual(values=c("black","grey")) + 
		ggplot2::xlab("y") + ggplot2::ylab("")

	return(p)
}

#' Plot Cross-Sections
#'
#' @export
#' @keywords internal
#' @param x a sstapreg object
#' @param stap_term name of stap term to plot
#' @param component one of c("Distance","Time")
#' @param fixed_val vector that contains fixed values for whichever component was not specified
#' @param p probability_interval
#' @param optional grid by default is NULL corresponding to a .1 step grid from the min-max of the covariate space
#'
plot_xsection <- function(x,stap_term = NULL, component = "Distance",fixed_val = 1, p = 0.95,grid = NULL)
	UseMethod("plot_xsection")

#' Plot Cross-Sections
#'
#' @export
#' @describeIn plot_xsection
#'
plot_xsection.sstapreg <- function(x,stap_term = NULL, component = "Distance",fixed_val =1 , p = 0.95,grid = NULL){

	Distance <- Time <- Median <- Grid <- Lower <- Upper <-  . <- .data <-  NULL
	check_p(p)
	spec <- x$specification
	if(!has_any_staps(spec))
		stop("Model has no stap terms specified")

	if(is.null(stap_term)){
		ix <- which(spec$component=="Distance-Time")[1]
		stap_term <- spec$term[ix]
	}
	else if(stap_term %in% spec$term){
		## check term is a stap_term
		stopifnot(get_component(spec,stap_term)=="Distance-Time")
	}else{
		stop("Term specified is not included in model")
	}

	beta <- as.matrix(x$stapfit)

	gd_eta <- get_stap(spec,stap_term,"Distance-Time",beta,x$family,grid)
	gd <- gd_eta$grid
	eta <- gd_eta$eta
	
	l <-  .5 - p/2
	u <- .5 + p/2
	
	ocomp <- switch(component,
	                "Distance"="Time",
	                "Time"="Distance")
	lbl <- paste0(ocomp, " fixed at ", fixed_val)
	
	ics <- which(gd[,ocomp]==fixed_val)
	gd <- gd[ics,]
	eta <- eta[ics,]

  	dplyr::tibble(Distance = gd$Distance,
				  Time = gd$Time,
  				  Lower = apply(eta,1,function(x) quantile(x,l)),
  				  Median = apply(eta,1,median),
  				  Upper = apply(eta,1,function(x) quantile(x,u))
				  ) -> pltdf


	if(has_bw(spec,stap_term)){

		dplyr::tibble(Distance = gd$Distance,
				Time = gd$Time,
				Lower = apply(gd_eta$eta_within[ics,],1,function(x) quantile(x,l)),
				Median = apply(gd_eta$eta_within[ics,],1,median),
				Upper = apply(gd_eta$eta_within[ics,],1,function(x) quantile(x,u)),
				Parameters = "within") -> pltdf2

		pltdf %>% dplyr::mutate(Parameters = "between") %>% rbind(.,pltdf2) -> pltdf

		pltdf %>% ggplot2::ggplot(ggplot2::aes(x=.data[[component]],y=Median)) + 
			ggplot2::geom_line()  + 
			ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),alpha=0.3) + 
			ggplot2::theme_bw() + 
			ggplot2::theme(strip.background=ggplot2::element_blank()) + 
			ggplot2::labs(subtitle = lbl) + 
			ggplot2::facet_wrap(~Parameters) -> pl2

		return(pl2)
	}

	pltdf %>% ggplot2::ggplot(ggplot2::aes(x = .data[[component]], y = Median )) + 
	ggplot2::geom_line() + 
	ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),alpha=0.3) + 
	ggplot2::theme_bw() + 
	ggplot2::geom_hline(ggplot2::aes(yintercept = 0),
						linetype=2,color='red') + 
	ggplot2::xlab(component) + 
	ggplot2::labs(subtitle = lbl) +
	ggplot2::ylab("") ->pl

	return(pl)

}

#' Plotting a random varying  coefficients regression model
#'
#' @export
#'
#' @param x rvcreg object
#'
plot.rvcreg <- function(x,...){

	grid <- dplyr::tibble(Distance = seq(from = floor(x$spec$ranges[1]),to = ceiling(x$spec$ranges[2]), by = 0.01))
	dmat <- mgcv::PredictMat(x$spec$smooth_objs[[1]],grid)
	beta <- get_coefnames(x$spec,'smooth')
	beta <- x$stanmat[,beta]
	eta <- tcrossprod(dmat,beta)
	pltdf <- get_pltdf(grid,eta,.95,"Distance")

	alpha <- get_coefnames(x$spec,)

	pltdf %>% 
		ggplot2::ggplot(ggplot2::aes(x = Distance,y = Median)) + 
	  ggplot2::geom_line() + 
		ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),alpha=0.3) + 
		ggplot2::theme_bw() + 
		ggplot2::geom_hline(ggplot2::aes(yintercept = 0),
							linetype=2,color='red') + 
		ggplot2::xlab("Distance") + 
		ggplot2::ylab("") ->pl

	return(pl)
}

# Internal ---------------------------------------------

check_p <- function(p){
	stopifnot(p<1 && p>0)
}

get_pltdf <- function(gd,eta,p,component,eta_wi = NULL){
	check_p(p)
	l <-  .5 - p/2
	u <- .5 + p/2

	pltdf <- dplyr::tibble(Lower = apply(eta,1,function(x) quantile(x,l)),
                         Median = apply(eta,1,median),
                         Upper = apply(eta, 1, function(x) quantile(x,u))
                         )

  pltdf <- cbind(pltdf,gd)

  if(!is.null(eta_wi)){
    pltdf <- pltdf %>% dplyr::mutate(Function = "Between")
    
    pltdf_  <- dplyr::tibble(Lower = apply(eta_wi,1, function(x) quantile(x,l)),
                  Median  = apply(eta_wi,1,median),
                  Function = "Within",
                  Upper = apply(eta_wi,1,function(x) quantile(x,u)))

    pltdf <- rbind(pltdf,
                  cbind(pltdf_,gd))

  }


	return(pltdf)

}


