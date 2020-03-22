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
#' @param formula Similar as for \code{\link[rstap]{stap_lm}}. 
#' @param subject_data required data argument
#' @param subject_id string containing the common id column in both data and distance data and/or time_data
#' @param basis_functions list with length equal to the number of BEFs in the stap_formula that specifies the basis function expansion
#' @param dt_data distance-time dataframe containing three columns (1) subj_ID, (2) BEF_name and (3) Distance between subj_ID and BEF 
#' @param BEF_col_name string name of the column in the dt dataframe that contains the character vectors describing the BEFs
#' @param distance_col_name string name of the column in the dt dataframe that contains the distance values
#' @param time_col_name string name of the column in the dt dataframe that contains the time values
#'
sstap_df <- function(stap_formula,
                     subject_data,
                     subject_id = NULL,
                     basis_functions = NULL,
                     dt_data = NULL,
                     BEF_col_name = NULL,
                     distance_col_name = NULL,
                     time_col_name = NULL){
  if(is.null(subject_id)){
    stop("subject_id must be provided so that distance and subject data can be joined")
  }
  if(!(BEF_col_name %in% colnames(dt_data)))
    stop("`BEF_col_name` must be a column in dt_data")
  if(!(distance_col_name %in% colnames(dt_data)) && !is.null(distance_col_name))
    stop("`distance_col_name` must be a column in dt_data")
  if(!(time_col_name %in% colnames(dt_data)) && !is.null(time_col_name))
    stop("`time_col_name` must be a column in dt_data")
  if(!(subject_id %in% colnames(dt_data)) || !(subject_id %in% colnames(subject_data)))
    stop("`subject_id` must be a column in both dt_data and subject_data")
  
  stap_data <- extract_stap_data(stap_formula)
  stcode <- stap_data$stap_code
  stlabels <- sapply(stcode,function(x) if(x%in%c(0,2)) return("Spatial") else if(x==1) return("Temporal"))
  group_labels <- sapply(seq_along(stap_data$group_indicator),function(x) if(stap_data$group_indicator[x]==1) return(paste0(stap_data$group_term[x],"_")) else return(""))
  BEFs <- stap_data$covariates
  subject_data <- subject_data %>% dplyr::arrange_(.dots=subject_id)
  dt_data <- dt_data %>% dplyr::arrange_(.dots=subject_id) %>% 
    dplyr::right_join(subject_data[,subject_id],by=subject_id)
  if(length(subject_id)>1)
    dt_data$new_id <- Reduce(function(x,y) stringr::str_c(x,y,sep="_"),purrr::map(subject_id,function(x) dt_data[,x,drop=TRUE]))
  else
    dt_data$new_id <- dt_data[,subject_id,drop=TRUE]
  dt_data$new_id <- factor(dt_data$new_id,levels=unique(dt_data$new_id))
  
  DirectEffect <- purrr::map(1:length(BEFs),function(ix){
    if(stcode[ix] %in% c(0,2))
      dt_col_name <- distance_col_name
    else
      dt_col_name <- time_col_name
    tmpdf <- dt_data %>% dplyr::filter(!!dplyr::sym(BEF_col_name)==BEFs[ix]|is.na(!!dplyr::sym(BEF_col_name)))
    df <- tmpdf %>% split(.[,"new_id"]) %>% 
      purrr::map(.,function(x){
        if(any(is.na(x[,dt_col_name,drop=TRUE])) ){
          return(c(0,basis_functions[[ix]](0)))
        }else if(stap_data$group_indicator[ix]==0){
          return(colSums(cbind(1,basis_functions[[ix]](x[,dt_col_name,drop=TRUE]))))
        }else{
          return(colSums(basis_functions[[ix]](x[,dt_col_name,drop=TRUE])))
        }
      })
    df <- do.call(rbind,df)
    colnames(df) <- paste0(stlabels[ix],"Effect_",group_labels[ix],BEFs[ix], "_",(0:(ncol(df)-1) + stap_data$group_indicator[ix]))
    return(df)
  })

  DirectEffect <- do.call(cbind,DirectEffect)
  if(any(stcode==2)){
    st_ics <- which(stcode==2)
    STimeDirectEffect <- purrr::map(st_ics,function(ix){
      tmpdf <- dt_data %>% dplyr::filter(!!dplyr::sym(BEF_col_name)==BEFs[ix])
      df <- tmpdf %>% split(.[,"new_id"]) %>% 
        purrr::map(.,function(x){
          dt_col_name <- time_col_name
          if(is.na(x[1,dt_col_name]))
            return(c(0,basis_functions[[ix]](0)))
          else 
            return(colSums(cbind(1,basis_functions[[ix]](x[,dt_col_name,drop=TRUE]))))
        })
      df <- do.call(rbind,df)
      colnames(df) <- paste0("TemporalEffect","_",group_labels[ix],BEFs[ix],(0:(ncol(df)-1) + stap_data$group_indicator[ix]))
      return(df)
    })
    STimeDirectEffect <- do.call(cbind,STimeDirectEffect)
    DirectEffect <- cbind(DirectEffect,STimeDirectEffect)
  }
  return(DirectEffect)
}
