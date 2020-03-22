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

#' Built Environment Networks Frequentist Generalized Linear Model 
#'
#'
#' @export
#'
#' @param formula Similar as for \code{\link[stats]{glm}}. 
#' @param stap_formula See \code{\link[rstap]{stap_glm}}
#' @param subject_data required data argument containing subject level outcome and covariates
#' @param subject_id string name for the common id column in both data and distance data and/or time_data
#' @param dt_data distance dataframe containing up to four columns:
#'  (1) subj_ID, (2) BEF_name and (3) Distance AND/OR (4) Time between subj_ID and BEF 
#' @param BEF_col_name string name for the column containing the BEF labels in dt_data dataframe
#' @param distance_col_name string name for the column containing the subject-BEF distances in the dt_data dataframe
#' @param time_col_name string name for the the column containing the subject-BEF times in the dt_data dataframe
#' @param family same as \code{\link[stats]{glm}}.
#' @param method one of c("glm","stan_glm" or "brm") for frequentist or bayesian method implimentation
#' @param ... args for \code{\link[stats]{glm}} or \code{\link[brms]{brm}} or \code{\link[rstanarm]{stan_glm}}
#' 
sstap_glm <- function(formula,
                    stap_formula,
                    subject_data,
                    subject_id = NULL,
                    basis_functions = NULL,
                    dt_data = NULL,
                    BEF_col_name = NULL,
                    distance_col_name = NULL,
                    time_col_name = NULL,
                    family = gaussian(),
					method = "brm",
                    ...){
	if(!(method %in% c('glm','stan_glm','brm')))
		stop("method must be one of c('glm','stan_glm','brm')")

	bef_df <- sstap_df(stap_formula = stap_formula,
					   subject_data = subject_data,
					   subject_id = subject_id,
					   basis_functions = basis_functions,
                       dt_data = dt_data,
                       BEF_col_name = BEF_col_name,
                       distance_col_name = distance_col_name,
                       time_col_name = time_col_name)
  
	if(family$family=="binomial"){
	if("cbind" %in% all.names(formula)){
	  resp <- paste("cbind(",all.vars(formula)[1],",",all.vars(formula)[2],")")
	  covs <- all.vars(formula)[3:length(all.vars(formula))]
	}
	}else{
	resp <- all.vars(formula)[1]
	covs <- all.vars(formula)[2:length(all.vars(formula))]
	}
  
  
	covs <- c(covs,colnames(bef_df))
	formula <- as.formula(paste(resp, " ~ ", paste(covs,collapse = " + ")))
	X <- cbind(subject_data,bef_df)
	if(method =="glm"){

		fit <- glm(formula,data = X, family = family,...)
		fit$basis_functions <- basis_functions
		fit$stap_data <- rstap:::extract_stap_data(stap_formula)
		fit$BEFs <- fit$stap_data$covariates
		if(any(fit$stap_data$stap_code %in% c(0,2)))
			fit$spaceranges <- lapply(fit$BEFs,function(x){ dt_data %>% 
				dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
				pull(!!dplyr::sym(distance_col_name)) %>% range(.)})
		if(any(fit$stap_data$stap_code %in% c(1,2)))
			fit$timeranges <- lapply(fit$BEFs,function(x){ dt_data %>% 
				dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
				pull(!!dplyr::sym(time_col_name)) %>% range(.)})

		structure(fit,class=c("glm","lm","sstap"))
	}else if(method == "stan_glm"){
		fit <- rstanarm::stan_glm(formula,data = X, family = family, ...)
		fit$basis_functions <- basis_functions
		fit$stap_data <- rstap:::extract_stap_data(stap_formula)
		fit$BEFs <- fit$stap_data$covariates
		if(any(fit$stap_data$stap_code %in% c(0,2)))
		fit$spaceranges <- lapply(fit$BEFs,function(x){ dt_data %>% 
			dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
			pull(!!dplyr::sym(distance_col_name)) %>% range(.)})
		if(any(fit$stap_data$stap_code %in% c(1,2)))
		fit$timeranges <- lapply(fit$BEFs,function(x){ dt_data %>% 
			dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
			pull(!!dplyr::sym(time_col_name)) %>% range(.)})
		fit$scales <- scales
		if(any(scales>1))
			warning("Note: BEF-covariates were scaled, see object$scale for scales used")
		structure(fit,class=c("stanreg","sstap"))
	}else if(method == "brms"){
		fit <- brms::brm(formula,data = X, family = family, ...)
		fit$basis_functions <- basis_functions
		fit$stap_data <- rstap:::extract_stap_data(stap_formula)
		fit$BEFs <- fit$stap_data$covariates
		if(any(fit$stap_data$stap_code %in% c(0,2)))
		fit$spaceranges <- lapply(fit$BEFs,function(x){ dt_data %>% 
			dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
			pull(!!dplyr::sym(distance_col_name)) %>% range(.)})
		if(any(fit$stap_data$stap_code %in% c(1,2)))
		fit$timeranges <- lapply(fit$BEFs,function(x){ dt_data %>% 
			dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
			pull(!!dplyr::sym(time_col_name)) %>% range(.)})
		fit$scales <- scales
		if(any(scales>1))
			warning("Note: BEF-covariates were scaled, see object$scale for scales used")

		structure(fit,class=c("brmsfit","sstap"))
	}
}
