#' @export
#' @rdname plot_effects
setGeneric("plot_effects_mer", function(object,pars = NULL) standardGeneric("plot_effects_mer"))

#' @export
#' @rdname plot_df 
setGeneric("plot_df_mer", function(object,pars = NULL) standardGeneric("plot_df_mer"))


#' Method for creating plot dataframe for bbMod objects
#' @export
#' 
setMethod("plot_df_mer",signature=signature("bbMod"),
	function(object,pars = NULL){
		if(!is.null(pars))
			code <- object@stap_code[which(object@BEFs %in% pars)]
		else
			code <- object@stap_code
		
		spatial_ics <- which(code == 0)
	  temporal_ics <- which(code == 1)
	  st_ics <- which(code == 2)
	  spatial_ics <- union(spatial_ics,st_ics)
	  temporal_ics <- union(temporal_ics,st_ics)
	  spatial_BEFs <- object@BEFs[spatial_ics]
	  temporal_BEFs <- object@BEFs[temporal_ics]
	  spatial_covs <- lapply(spatial_BEFs,function(x) grep(paste0(x,"_"),names(lme4::fixef(object)),value=TRUE))
	  temporal_covs <- lapply(temporal_BEFs,function(x) grep(paste0(x,"_"),names(lme4::fixef(object)),value=TRUE))
	  spacegrids <- lapply(spatial_ics,function(x) seq(from = object@spaceranges[[x]][1],
	                                                   to = object@spaceranges[[x]][2],
	                                                   by = 0.01))
	  timegrids <- lapply(temporal_ics,function(x) seq(from = object@timeranges[[x]][1],
	                                                   to = object@timeranges[[x]][2],
	                                                   by = 0.01))
	  spacegridmats <- purrr::map(seq_along(spacegrids),function(x) cbind(1,predict(object@basis_functions[[x]](0),spacegrids[[x]])))
	  timegridmats <- purrr::map(seq_along(timegrids),function(x) cbind(1,predict(object@basis_functions[[x]](0),timegrids[[x]])))
      pltdf <- pltdf_helper(object,spacegrids,spacegridmats,spatial_covs,spatial_BEFs,
	                      timegridmats,temporal_covs,temporal_BEFs)
      return(pltdf)
})

#' Method for plotting bbMod objects
#' @export
#' 
setMethod("plot_effects_mer",signature = signature("bbMod"),
			function(object,pars= NULL){
			  pltdf <- plot_df_mer(object,pars)
			  p <- pltdf %>% ggplot2::ggplot(ggplot2::aes(x=Grid,y=Effect)) + 
				ggplot2::geom_line() +
				ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),alpha=0.3) + 
				ggplot2::facet_wrap(~label) + 
				ggplot2::theme(strip.background=ggplot2::element_blank()) +  
				ggplot2::labs(title = "Exposure Effect",
							  subtitle = "Shaded Area represents 95% Confidence Interval")
				return(p)
			  }
)
