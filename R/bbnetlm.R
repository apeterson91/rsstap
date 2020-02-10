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

#' Built Environment Networks Linear Model 
#'
#'
#' @export
#'
#' @param formula Similar as for \code{\link[stats]{lm}}. 
#' @param data required data argument
#' @param subject_id string containing the common id column in both data and distance data and/or time_data
#' @param BEFs character vector containing BEF's coded in distance/time dataframe to include
#' @param distance_data distance dataframe containing up to four columns:
#'  (1) subj_ID, (2) BEF_name and (3) Distance OR (4) Time between subj_ID and BEF 
#' @param ... args for \code{\link[stats]{lm}}
#' 
bbnet_lm <- function(formula,
                     subject_data,
                     subject_id = NULL,
                     BEFs = NULL,
                     basis_functions = NULL,
                     dt_data = NULL,
                     BEF_col_name = NULL,
                     dt_col_name = NULL,
                     ...){
  
  if(is.null(BEFs))
    BEFs <- unique(dt_data[,BEF_col_name])
  
  bef_df <- bbnet_df(subject_data,
                     subject_id,
                     BEFs,
                     basis_functions,
                     dt_data,
                     BEF_col_name,
                     dt_col_name)
  scales <- apply(bef_df,2,function(x) round(median(x),2))
  bef_df <- bef_df / scales
  resp <- all.vars(formula)[1]
  covs <- all.vars(formula)[2:length(all.vars(formula))]
  covs <- c(covs,colnames(bef_df))
  formula <- as.formula(paste(resp, " ~ ", paste(covs,collapse = " + ")))
  X <- cbind(subject_data,bef_df)
  
  fit <- lm(formula,data = X, ...)
  
  fit$coefficients[colnames(bef_df)] <- fit$coef[colnames(bef_df)] / scales
  fit$basis_functions <- basis_functions
  fit$ranges <- lapply(BEFs,function(x){ dt_data %>% 
                         dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
                         pull(!!dplyr::sym(dt_col_name)) %>% range(.)})
  fit$BEFs <- BEFs
  
  structure(fit,class=c("lm","bbnet"))
}
