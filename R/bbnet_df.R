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

#' Built Environment Network Direct Effect Matrix
#'
#'
#' @export
#'
#' @param formula Similar as for \code{\link[stats]{lm}}. 
#' @param subject_data required data argument
#' @param subject_id string containing the common id column in both data and distance data and/or time_data
#' @param dt_data distance-time dataframe containing three columns (1) subj_ID, (2) BEF_name and (3) Distance between subj_ID and BEF 
#' @param BEF_col_name string name of the column in the dt dataframe that contains the character vectors describing the BEFs
#' @param distance_col_name string name of the column in the dt dataframe that contains the distance values
#' @param time_col_name string name of the column in the dt dataframe that contains the time values
#'
bbnet_df <- function(stap_formula,
                     subject_data,
                     subject_id = NULL,
                     basis_functions = NULL,
                     dt_data = NULL,
                     BEF_col_name = NULL,
                     distance_col_name = NULL,
                     time_col_name = NULL){
  stap_data <- rstap:::extract_stap_data(stap_formula)
  stcode <- stap_data$stap_code
  stlabels <- sapply(stcode,function(x) if(x%in%c(0,2))return("Spatial") else if(x==1) return("Temporal"))
  BEFs <- stap_data$covariates
  subject_data <- subject_data %>% dplyr::arrange(!!dplyr::sym(subject_id))
  dt_data <- dt_data %>% dplyr::arrange(!!dplyr::sym(subject_id)) %>% 
    right_join(subject_data[,subject_id],by=subject_id)
  DirectEffect <- purrr::map(1:length(BEFs),function(ix){
    tmpdf <- dt_data %>% dplyr::filter(!!dplyr::sym(BEF_col_name)==BEFs[ix])
    df <- tmpdf %>% split(.[,subject_id]) %>% 
      purrr::map(.,function(x){
        if(stcode[ix] %in% c(0,2))
          dt_col_name <- distance_col_name
        else
          dt_col_name <- time_col_name
        if(is.na(x[1,dt_col_name]))
          return(c(0,basis_functions[[ix]](0)))
        else
          return(colSums(cbind(1,basis_functions[[ix]](x[,dt_col_name,drop=TRUE]))))
      })
    df <- do.call(rbind,df)
    colnames(df) <- paste0("Direct",stlabels[ix],"Effect_",BEFs[ix], "_",0:(ncol(df)-1))
    return(df)
  })
  DirectEffect <- do.call(cbind,DirectEffect)
  if(any(stcode==2)){
    st_ics <- which(stcode==2)
    STimeDirectEffect <- purrr::map(st_ics,function(ix){
      tmpdf <- dt_data %>% dplyr::filter(!!dplyr::sym(BEF_col_name)==BEFs[ix])
      df <- tmpdf %>% split(.[,subject_id]) %>% 
        purrr::map(.,function(x){
          dt_col_name <- time_col_name
          if(is.na(x[1,dt_col_name]))
            return(c(0,basis_functions[[ix]](0)))
          else
            return(colSums(cbind(1,basis_functions[[ix]](x[,dt_col_name,drop=TRUE]))))
        })
      df <- do.call(rbind,df)
      colnames(df) <- paste0("DirectTemporalEffect","_",0:(ncol(ST)))
      return(df)
    })
    STimeDirectEffect <- do.call(cbind,STimeDirectEffect)
    DirectEffect <- cbind(DirectEffect,STimeDirectEffect)
  }
  return(DirectEffect)
}
