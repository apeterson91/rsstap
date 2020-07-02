

#' Plot Spline Spatial Temporal Aggregated Predictors
#'
#' @export
#' @param x sstapreg object
#' 
plot.sstapreg <- function(x,stap_term = NULL,component = NULL){

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
	if(x$stap_components[1] == "Distance") 
		gd <- data.frame(Distance = seq(from = floor(x$ranges[[ix]]$Distance[1]), to = ceiling(x$ranges[[ix]]$Distance[2]), by =0.01))
	else if(x$stap_components[1] == "Time") 
		gd <- data.frame(Time = seq(from = floor(x$ranges[[ix]]$Time[1]), to = ceiling(x$ranges[[ix]]$Time[2]), by =0.01))
	else if(x$stan_components[1] == "Distance-Time")
		gd <- data.frame(Distance = seq(from = floor(bef$ranges[[ix]]$Distance[1]), to = ceiling(bef$ranges[[ix]]$Distance[2]), by =0.01),
						 Time = seq(from = floor(bef$ranges[[ix]]$Time[1]), to = ceiling(bef$ranges[[ix]]$Time[2]), by =0.01))

	mat <- mgcv::Predict.matrix(bef,gd)

	beta <- as.matrix(x$stapfit)
	betas <- grep("beta",colnames(beta))
	beta <- beta[,betas]
	beta <- beta[,x$ind[[ix]]]
	eta <- tcrossprod(mat,beta)
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
		ggplot2::geom_rug(sides='b') + 
		ggplot2::ylab("")
	
}

#' Posterior Predictive Checks 
#'
#' @export
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
						  tidyr::gather(contains('yhat'),key="Parameter",value="Samples") %>% 
						  dplyr::mutate(Parameter = 'yrep')
	pltdf <- rbind(pltdf,
				   dplyr::tibble(iteration_ix = 0, Parameter='y',Samples= x$model$y ))
	
	p <- pltdf %>% 
		ggplot2::ggplot(ggplot2::aes(x=Samples,color=Parameter,group=iteration_ix)) + 
		ggplot2::geom_density() + ggplot2::theme_bw()+ ggplot2::theme(legend.title=ggplot2::element_blank()) +   
		ggplot2::scale_colour_manual(values=c("black","grey")) + 
		xlab("y") + ylab("")

	return(p)
}
