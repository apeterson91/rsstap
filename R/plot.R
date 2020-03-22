# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3 # of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Parent Method for plot_directeffect
#' @export
#' 
plot_effects <- function(object,pars = NULL){
  UseMethod("plot_effects")
}

#' Parent Method for getting plot dataframe
#' @export
#' 
plot_df <- function(object,pars = NULL){
  UseMethod("plot_df")
}

#' Plots the Direct Spatio-Temporal Exposure of a sstap Model
#'
#'
#' @export
#'
#' @param object of type sstap
#' @param pars optional parameter selection
#' 
plot_effects.sstap <- function(object,pars = NULL){
  
  pltdf <- plot_df(object,pars)
  if("lm" %in% class(object) || "lmerMod" %in% class(object))
    subtitle_string <- "Shaded Area represents 95% Confidence Interval"
  else
    subtitle_string <- "Shaded Area represents 95% Credible Interval"
  p <- pltdf %>% ggplot2::ggplot(ggplot2::aes(x=Grid,y=Effect)) + 
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),alpha=0.3) + 
    ggplot2::facet_wrap(~label) + ggplot2::theme(strip.background=ggplot2::element_blank()) +  
    ggplot2::labs(title = "Exposure Effect",
                  subtitle = subtitle_string)
    return(p)
}

#'
#'@export  
plot_df.sstap <- function(object,pars){
  
  if(!is.null(pars))
    code <- object$stap_data$stap_code[which(object$stap_data$covariates %in% pars)]
  else 
    code <- object$stap_data$stap_code
  spatial_ics <- which(code==0)
  temporal_ics <- which(code==1)
  st_ics <- which(code==2)
  spatial_ics <- union(spatial_ics,st_ics)
  temporal_ics <- union(temporal_ics,st_ics)
  spatial_BEFs <- object$BEFs[spatial_ics]
  temporal_BEFs <- object$BEFs[temporal_ics]
  spatial_BEFs_grp <- get_grp_string(spatial_BEFs)
  temporal_BEFs_grp <- get_grp_string(temporal_BEFs)
  if(length(intersect(c("lm","stanreg"),class(object)))>0){
    spatial_covs <- lapply(spatial_BEFs_grp,function(x) grep(x,names(coef(object)),value=TRUE))
    temporal_covs <- lapply(temporal_BEFs_grp,function(x) grep(x,names(coef(object)),value=TRUE))  
  }
  else if("brmsfit" %in% class(object)){
    spatial_covs <- lapply(spatial_BEFs_grp,function(x) paste0("b_",grep(x,rownames(brms::fixef(object)),value=TRUE)))
    temporal_covs <- lapply(temporal_BEFs_grp,function(x) paste0("b_",grep(x,rownames(brms::fixef(object)),value=TRUE)))
  }
  spacegrids <- lapply(spatial_ics,function(x) seq(from = object$spaceranges[[x]][1],
                                                   to = object$spaceranges[[x]][2],
                                                   by = 0.01))
  timegrids <- lapply(temporal_ics,function(x) seq(from = object$timeranges[[x]][1],
                                                   to = object$timeranges[[x]][2],
                                                   by = 0.01))
  spacegridmats <- purrr::map(seq_along(spacegrids),function(x) cbind(1,predict(object$basis_functions[[x]](0),spacegrids[[x]])))
  timegridmats <- purrr::map(seq_along(timegrids),function(x) cbind(1,predict(object$basis_functions[[x]](0),timegrids[[x]])))
	pltdf <- pltdf_helper(object,spacegrids,spacegridmats,spatial_covs,spatial_BEFs,
	                      timegridmats,temporal_covs,temporal_BEFs)
	
  return(pltdf)
}

pltdf_helper <- function(object,spacegrids,spacegridmats,spatial_covs,spatial_BEFs,
                         timegrids,timegridmats,temporal_covs,temporal_BEFs){

	if("lm" %in% class(object) || "sstapMod" %in% class(object)){
	  intervals <-  confint(object)
	  if("sstapMod" %in% class(object))
	    coefs <- lme4::fixef(object)
	  else
	    coefs <- coef(object)
	}
	else if("stanreg" %in% class(object) || "brmsfit" %in% class(object)){
	  intervals <- rstanarm::posterior_interval(object)
	  if("brmsfit" %in% class(object)){
	    coefs <- brms::fixef(object)[,1]
	    names(coefs) <- paste0("b_",names(coefs))
	  }
	  else
	    coefs <- coef(object)
	}
    effects <- lapply(seq_along(spacegridmats),function(x) spacegridmats[[x]] %*% coefs[spatial_covs[[x]]])
    lowers <- lapply(seq_along(spacegridmats),function(x) spacegridmats[[x]] %*% intervals[spatial_covs[[x]],1])
    uppers <- lapply(seq_along(spacegridmats),function(x) spacegridmats[[x]] %*% intervals[spatial_covs[[x]],2])
    pltdf <- purrr::map_dfr(seq_along(effects),function(x) dplyr::tibble(Grid = spacegrids[[x]],
                                                           Effect = effects[[x]],
                                                           lower = lowers[[x]],
                                                           upper = uppers[[x]],
                                                           BEF = spatial_BEFs[x],
                                                           type = "Spatial",
                                                           label = paste(type,BEF)))
    effects <- lapply(seq_along(timegridmats),function(x) timegridmats[[x]] %*% coefs[temporal_covs[[x]]])
    lowers <- lapply(seq_along(timegridmats),function(x) timegridmats[[x]] %*% intervals[temporal_covs[[x]],1])
    uppers <- lapply(seq_along(timegridmats),function(x) timegridmats[[x]] %*% intervals[temporal_covs[[x]],2])
    pltdf <- rbind(pltdf,purrr::map_dfr(seq_along(effects),function(x) dplyr::tibble(Grid = timegrids[[x]],
                                                                        Effect = effects[[x]],
                                                                        lower = lowers[[x]],
                                                                        upper = uppers[[x]],
                                                                        BEF = temporal_BEFs[x],
                                                                        type = "Temporal",
                                                                        label = paste(type,BEF))))

	return(pltdf)
}

#' pastes string for regex expression
get_grp_string <- function(BEFs){
  if(length(BEFs)==0)
    return(character())
  else
    return(paste0(BEFs,"_[0-9]"))
}
