

#' Posterior Predictive Distribution Estimate
#'
#' @export
#' 
#' @aliases posterior_predict
#' @param object sstapreg object
#' @param type one of "link","response". See \code{\link[stats]{predict.glm}}.
#' @param ... currently ignored
#' 
posterior_predict.sstapreg <- function(object,type = 'link',...){

	nms <- c(colnames(object$model$Z),Reduce(union,lapply(object$specification$X,function(y) colnames(y))))
	pars <- as.matrix(object)[,nms]
	mat <- cbind(object$model$Z,Reduce(cbind,object$specification$X))
	eta <- t(tcrossprod(mat,pars))
	if(type %in% c("response"))
		eta <- object$family$linkinv(eta)
	return(eta)
}

