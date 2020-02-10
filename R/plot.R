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
plot_directeffect <- function(object,pars = NULL){
  UseMethod("plot_directeffect")
}

#' Parent Method for getting plot dataframe
#' @export
#' 
plot_df <- function(object,pars = NULL){
  UseMethod("plot_df")
}

#' Plots the Direct Spatio-Temporal Exposure of a bbnet Model
#'
#'
#' @export
#'
#' @param object of type bbnet
#' @param pars optional parameter selection
#' 
plot_directeffect.bbnet <- function(object,pars = NULL){
  
  pltdf <- plot_df(object,pars)
  if("lm" %in% class(object))
    subtitle_string <- "Shaded Area represents 95% Asymptotic Confidence Interval"
  else
    subtitle_string <- "Shaded Area represents 90% Credible Interval"
  p <- pltdf %>% ggplot2::ggplot(ggplot2::aes(x=Grid,y=Effect)) + 
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),alpha=0.3) + 
    ggplot2::facet_wrap(~BEF) + 
    ggplot2::labs(title = "Exposure Effect",
                  subtitle = subtitle_string)
    return(p)
}

#'
#'@export 
plot_df.bbnet <- function(object,pars){
  
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
  spatial_covs <- lapply(spatial_BEFs,function(x) grep(x,names(coef(object)),value=TRUE))
  temporal_BEFs <- object$BEFs[temporal_ics]
  temporal_covs <- lapply(temporal_BEFs,function(x) grep(x,names(coef(object)),value=TRUE))
  spacegrids <- lapply(spatial_ics,function(x) seq(from = object$spaceranges[[x]][1],
                                                   to = object$spaceranges[[x]][2],
                                                   by = 0.01))
  timegrids <- lapply(temporal_ics,function(x) seq(from = object$timeranges[[x]][1],
                                                   to = object$timeranges[[x]][2],
                                                   by = 0.01))
  spacegridmats <- purrr::map(seq_along(spacegrids),function(x) cbind(1,predict(object$basis_functions[[x]](0),spacegrids[[x]])))
  timegridmats <- purrr::map(seq_along(timegrids),function(x) cbind(1,predict(object$basis_functions[[x]](0),timegrids[[x]])))
  if("lm" %in% class(object)){
    effects <- lapply(seq_along(spacegridmats),function(x) spacegridmats[[x]] %*% coef(object)[spatial_covs[[x]]])
    lowers <- lapply(seq_along(spacegridmats),function(x) spacegridmats[[x]] %*% confint(object)[spatial_covs[[x]],1])
    uppers <- lapply(seq_along(spacegridmats),function(x) spacegridmats[[x]] %*% confint(object)[spatial_covs[[x]],2])
    pltdf <- purrr::map_dfr(seq_along(effects),function(x) tibble(Grid = spacegrids[[x]],
                                                           Effect = effects[[x]],
                                                           lower = lowers[[x]],
                                                           upper = uppers[[x]],
                                                           BEF = spatial_BEFs[x],
                                                           type = "Spatial",
                                                           label = "Space"))
    effects <- lapply(seq_along(timegridmats),function(x) timegridmats[[x]] %*% coef(object)[temporal_covs[[x]]])
    lowers <- lapply(seq_along(timegridmats),function(x) timegridmats[[x]] %*% confint(object)[temporal_covs[[x]],1])
    uppers <- lapply(seq_along(timegridmats),function(x) timegridmats[[x]] %*% confint(object)[temporal_covs[[x]],2])
    pltdf <- rbind(pltdf,purrr::map_dfr(seq_along(effects),function(x) tibble(Grid = timegrids[[x]],
                                                                        Effect = effects[[x]],
                                                                        lower = lowers[[x]],
                                                                        upper = uppers[[x]],
                                                                        BEF = temporal_covs[x],
                                                                        type = "Temporal",
                                                                        label = "Time")))
  }else if("stanreg" %in% class(object)){
    effects <- lapply(seq_along(spacegridmats),function(x) spacegridmats[[x]] %*% coef(object)[spatial_covs[[x]]])
    lowers <- lapply(seq_along(spacegridmats),function(x) spacegridmats[[x]] %*% rstanarm::posterior_interval(object)[spatial_covs[[x]],1])
    uppers <- lapply(seq_along(spacegridmats),function(x) spacegridmats[[x]] %*% rstanarm::posterior_interval(object)[spatial_covs[[x]],2])
    pltdf <- purrr::map_dfr(seq_along(effects),function(x) tibble(Grid = spacegrids[[x]],
                                                   Effect = effects[[x]],
                                                   lower = lowers[[x]],
                                                   upper = uppers[[x]],
                                                   BEF = object$BEFs[x],
                                                   type = "Spatial",
                                                   label = "Distance"))
    effects <- lapply(seq_along(timegridmats),function(x) timegridmats[[x]] %*% coef(object)[temporal_covs[[x]]])
    lowers <- lapply(seq_along(timegridmats),function(x) timegridmats[[x]] %*% rstanarm::posterior_interval(object)[temporal_covs[[x]],1])
    uppers <- lapply(seq_along(timegridmats),function(x) timegridmats[[x]] %*% rstanarm::posterior_interval(object)[temporal_covs[[x]],2])
    pltdf <- rbind(pltdf,purrr::map_dfr(temporal_ics,function(x) tibble(Grid = timegrids[[x]],
                                                                        Effect = effects[[x]],
                                                                        lower = lowers[[x]],
                                                                        upper = uppers[[x]],
                                                                        BEF = temporal_covs[x],
                                                                        type = "Temporal",
                                                                        label = "Time")))
  }
  return(pltdf)
}
