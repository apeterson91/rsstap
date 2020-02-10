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

#' Parent Method for plot_effect
#' @export
plot_effect <- function(object,pars = NULL){
  UseMethod("plot_effect")
}

#' Plotting Methods for bbnet objects
#'
#'
#' @export
#'
#' @param object of type bbnet
#' @param pars optional parameter selection
plot_effect.bbnet <- function(object,pars = NULL){
  
  q <- length(object$BEFs)
  grids <- lapply(object$spaceranges,function(x) seq(from = x[1],to=x[2],by=0.01))
  grid_mats <- purrr::map(1:length(grids),function(x) cbind(1,predict(object$basis_functions[[x]](0),grids[[x]])))
  coef_names <- lapply(object$BEFs,function(x) grep(x,names(coef(object)),value=TRUE))
  if("lm" %in% class(object)){
    effects <- lapply(1:q,function(x) grid_mats[[x]] %*% coef(object)[coef_names[[x]]])
    lowers <- lapply(1:q,function(x) grid_mats[[x]] %*% confint(object)[coef_names[[x]],1])
    uppers <- lapply(1:q,function(x) grid_mats[[x]] %*% confint(object)[coef_names[[x]],2])
    pltdf <- purrr::map_dfr(1:q,function(x) tibble(Distance = grids[[x]],
                                                   Effect = effects[[x]],
                                                   lower = lowers[[x]],
                                                   upper = uppers[[x]],
                                                   BEF = object$BEFs[x]))
    p <- pltdf %>% ggplot2::ggplot(ggplot2::aes(x=Distance,y=Effect)) + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),alpha=0.3) +
      ggplot2::geom_line() + ggplot2::facet_wrap(~BEF)
    return(p)
  }else if("stanreg" %in% class(object)){
    effects <- lapply(1:q,function(x) grid_mats[[x]] %*% coef(object)[coef_names[[x]]])
    lowers <- lapply(1:q,function(x) grid_mats[[x]] %*% rstanarm::posterior_interval(object)[coef_names[[x]],1])
    uppers <- lapply(1:q,function(x) grid_mats[[x]] %*% rstanarm::posterior_interval(object)[coef_names[[x]],2])
    pltdf <- purrr::map_dfr(1:q,function(x) tibble(Distance = grids[[x]],
                                                   Effect = effects[[x]],
                                                   lower = lowers[[x]],
                                                   upper = uppers[[x]],
                                                   BEF = object$BEFs[x]))
    p <- pltdf %>% ggplot2::ggplot(ggplot2::aes(x=Distance,y=Effect)) + 
      ggplot2::geom_line() + 
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower,ymax=upper),
                                                  alpha=0.3) + 
      ggplot2::facet_wrap(~BEF)
    return(p)
    
    
  }
  

}
