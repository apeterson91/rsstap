

#' Plot Spline Spatial Temporal Aggregated Predictors
#'
#' @export
#' @param x sstapreg object
#' @param stap_term optional string for name of BEF smooth function to plot. 
#' Alternatively plots first BEF smooth function 
#' @param component one of c("Distance","Time","Distance-Time") 
#' corresponding to the smooth function domain
#' @param ... ignored
#' 
plot.sstapreg <- function(x,stap_term = NULL,component = NULL,...){

		Distance <- Time <- Median <- Grid <- Lower <- Upper <-  NULL
	if(is.null(stap_term)){
		ix <- 1
		stap_term <- x$stap_terms[ix]
		component <- x$stap_components[ix]
		stopifnot(x$stap_components[ix] == component)
		stopifnot(x$stap_term[ix] == stap_term)
	}
	else if(stap_term %in% x$stap_terms){
		ix <- which(stap_term == x$stap_terms)
		ix <- intersect(ix,which(x$stap_components==component))
	}
	else
		stop("stap_term must be NULL or one of the stap terms in the model object")
	bef <- x$model$smooths[[ix]][[1]]
	if(component == "Distance") 
		gd <- data.frame(Distance = seq(from = floor(x$ranges[[ix]]$Distance[1]), to = ceiling(x$ranges[[ix]]$Distance[2]), by =0.01))
	else if(component == "Time") 
		gd <- data.frame(Time = seq(from = floor(x$ranges[[ix]]$Time[1]), to = ceiling(x$ranges[[ix]]$Time[2]), by =0.01))
	else if(component == "Distance-Time"){
		gd <- expand.grid(Distance = seq(from = floor(x$ranges[[ix]]$Distance[1]), to = ceiling(x$ranges[[ix]]$Distance[2]), by=0.01),
						  Time = seq(from = floor(x$ranges[[ix]]$Time[1]), to = ceiling(x$ranges[[ix]]$Time[2]), by=0.01))
		gd <- as.data.frame(gd)
	}
		

	mat <- mgcv::Predict.matrix(bef,gd)

	beta <- as.matrix(x$stapfit)
	betas <- grep("beta",colnames(beta))
	beta <- beta[,betas]
	beta <- beta[,x$ind[[ix]]]
	eta <- tcrossprod(mat,beta)
	if(component %in% c("Distance","Time")){
  	dplyr::tibble(Grid = gd[,1],
  				  Lower = apply(eta,1,function(x) quantile(x,0.025)),
  				  Median = apply(eta,1,median),
  				  Upper = apply(eta,1,function(x) quantile(x,0.975))) -> pltdf
  
  	pltdf %>% ggplot2::ggplot(ggplot2::aes(x = Grid, y = Median )) + 
  		ggplot2::geom_line() + 
  		ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),alpha=0.3) + 
  		ggplot2::theme_bw() + ggplot2::geom_hline(ggplot2::aes(yintercept = 0),
  												  linetype=2,color='red') + 
  		ggplot2::xlab(component) + 
  		ggplot2::ylab("") ->pl
	}else{
	  dplyr::tibble(Distance = gd$Distance,
	                Time = gd$Time,
	                Median = apply(eta,1,median)) -> pltdf
	  pltdf %>% ggplot2::ggplot(ggplot2::aes(x=Distance,y=Time,z=Median)) + 
	    ggplot2::geom_contour()  ->pl
	}
	
	return(pl)
}


#' 3D plots for rsstap models
#'
#' @export
#' @param x sstapreg object
#' @param stap_term string argument for which bef term to plot
#'
plot3D <- function(x,stap_term = NULL)
	UseMethod("plot3D")

#'
#' @describeIn plot3D 3d stap plot
#' @export 
#'
plot3D.sstapreg <- function(x,stap_term = NULL){

	if(is.null(stap_term)){
		ix <- which(x$stap_components == "Distance-Time")
		if(length(ix)<=0)
			stop("No Tensor-Spatial Temporal Terms included in Model")
		ix <- ix[1]
		stap_term <- x$stap_terms[ix]
		component <- x$stap_components[ix]
		stopifnot(x$stap_components[ix] == component)
		stopifnot(x$stap_term[ix] == stap_term)
	}
	else if(stap_term %in% x$stap_terms){
		ix <- which(stap_term == x$stap_terms)
		ix <- intersect(ix,which(x$stap_components==component))
		if(length(ix)<=0)
			stop("No Tensor-Spatial Temporal Terms included in Model")
	}
	bef <- x$model$smooths[[ix]][[1]]

	gd <- expand.grid(Distance = seq(from = floor(x$ranges[[ix]]$Distance[1]), 
									 to = ceiling(x$ranges[[ix]]$Distance[2]), by=0.01),
					  Time = seq(from = floor(x$ranges[[ix]]$Time[1]),
								 to = ceiling(x$ranges[[ix]]$Time[2]), by=0.01))
	gd <- as.data.frame(gd)
		

	mat <- mgcv::Predict.matrix(bef,gd)

	beta <- as.matrix(x$stapfit)
	betas <- grep("beta",colnames(beta))
	beta <- beta[,betas]
	beta <- beta[,x$ind[[ix]]]
	eta <- tcrossprod(mat,beta)
  	dplyr::tibble(Distance = gd$Distance,
				  Time = gd$Time,
  				  Lower = apply(eta,1,function(x) quantile(x,0.025)),
  				  Exposure = apply(eta,1,median),
  				  Upper = apply(eta,1,function(x) quantile(x,0.975))) -> pltdf
  
	pl <- plotly::plot_ly(pltdf,x = ~Distance, y = ~Time, z = ~Exposure,
						  mode = "markers", 
						  type='scatter3d') 

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

	samp <- sample(1:nrow(as.matrix(x$stapfit)),num_reps)
	yhatmat <- as.matrix(x$stapfit)
	yhats <- grep("yhat",colnames(yhatmat))
	yhatmat <- yhatmat[samp,yhats]
	pltdf <- suppressMessages(dplyr::as_tibble(yhatmat)) %>% dplyr::mutate(iteration_ix = 1:dplyr::n()) %>%
						  tidyr::gather(dplyr::contains('yhat'),key="Parameter",value="Samples") %>% 
						  dplyr::mutate(Parameter = 'yrep')
	if(is.matrix(x$model$y))
	  y <- x$model$y[,1]
	else
	  y <- x$model$y
	
	pltdf <- rbind(pltdf,
				   dplyr::tibble(iteration_ix = 0, Parameter='y',Samples= y ))
	
	p <- pltdf %>% 
		ggplot2::ggplot(ggplot2::aes(x=Samples,color=Parameter,group=iteration_ix)) + 
		ggplot2::geom_density() + ggplot2::theme_bw()+ ggplot2::theme(legend.title=ggplot2::element_blank()) +   
		ggplot2::scale_colour_manual(values=c("black","grey")) + 
		ggplot2::xlab("y") + ggplot2::ylab("")

	return(p)
}
