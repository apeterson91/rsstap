#' Randomly Varying Coefficients Spatial Temporal Aggregated Predictor Model
#'
#' Fits a spatial temporal aggregated predictor model where the 
#' spatial-temporal function is modeled using splines and
#' the function is modeled non additively as a function of proximate BEFs
#' 
#' @export
#'
#' @param mean_formula formula for mean function
#' @param rvc_formula formula for random varying coefficients
#' @param benvo built environment object from the rbenvo package containing the relevant data
#' @param weights for unequal variances
#' @param ... arguments for stan sampler
rvc_lm <- function(mean_formula,
				   rvc_formula,
				   benvo,
				   QR = TRUE,
				   weights = NULL,
				   ...){

	spec <- get_rvcspec(mean_formula,
						rvc_formula,
						benvo)

	sparts <- rstan::extract_sparse_parts(spec$AggMat)
	
	standata <- list(N = length(spec$mf$y),
					 P = ncol(spec$mf$X),
					 L = ncol(spec$X),
					 M = nrow(spec$X),
					 K = ncol(spec$D),
					 w_num = length(sparts$w), 
					 v_num = length(sparts$v),
					 v = sparts$v, 
					 w = sparts$w,
					 u = sparts$u,
					 y = spec$mf$y,
					 X_smooth = spec$X,
					 X_fixef = spec$mf$X,
					 S1 = spec$S1,
					 S2 = spec$S2,
					 D = spec$D,
					 has_weights = 0,
					 weights = as.numeric()
	)
	if(!is.null(weights)){
		standata$weights <- weights
		standata$has_weights <- 1
	}

	sampling_args <- set_sampling_args(
						object = stanmodels$rvc_continuous,
						control = list(adapt_delta = 0.85,
						               max_treedepth = 10),
						pars = c("beta_fixef","beta_smooth","alpha","sigma","tau1","tau2"),
						data = standata,
						show_messages = FALSE,
						save_warmup = FALSE,
						...
						)

	fit <- do.call(sampling,sampling_args)

	new_names <- c(get_coefnames(spec),
				   "sigma",
					paste0("smooth_precision[tau",1:2,"]"),
				   "log-posterior"
					)

    fit@sim$fnames_oi <- new_names

	out <- list(rvcfit = fit,
				spec = spec,
				mean_formula = mean_formula,
				rvc_formula = rvc_formula
			)

	return(rvcreg(out))


}
