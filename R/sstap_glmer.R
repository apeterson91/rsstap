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

#' Built Environment Networks Generalized Linear Mixed Effects Regression Model 
#'
#'
#' @export
#'
#' @param formula Similar as for \code{\link[lme4]{lmer}}. 
#' @param stap_formula See \code{\link[rstap]{stap_lm}}
#' @param subject_data required data argument containing subject level outcome and covariates
#' @param subject_id string name for the common id column in both data and distance data and/or time_data
#' @param basis_functions  list of basis functions for modeling the spatial (and/or) temporal exposure to BEFs
#' @param dt_data distance dataframe containing up to four columns:
#'  (1) subj_ID, (2) BEF_name and (3) Distance AND/OR (4) Time between subj_ID and BEF 
#' @param BEF_col_name string name for the column containing the BEF labels in dt_data dataframe
#' @param distance_col_name string name for the column containing the subject-BEF distances in the dt_data dataframe
#' @param time_col_name string name for the the column containing the subject-BEF times in the dt_data dataframe
#' @param family same as \code{\link[stats]{glm}}.
#' @param method one of c("glmer","stan_glmer" or "brm") for frequentist or bayesian method implimentation
#' @param ... args for \code{\link[stats]{lm}}
#' @importFrom methods new
#' 
sstap_glmer <- function(formula,
                    stap_formula,
                    subject_data,
                    subject_id = NULL,
                    basis_functions = NULL,
                    dt_data = NULL,
                    BEF_col_name = NULL,
                    distance_col_name = NULL,
                    time_col_name = NULL,
					family = stats::gaussian(),
					method = "brm",
                    ...){
	if(!(method %in% c('glmer','stan_glmer','brm')))
		stop("method must be one of c('glmer','stan_glmer','brm')")
	if(any(is.null(c(subject_id,basis_functions,dt_data,BEF_col_name))))
		stop("subject_id, basis_functions, dt_data and BEF_col_name cannot be NULL")
	if(is.null(distance_col_name) && is.null(time_col_name))
		stop("At least one of distance_col_name or time_col_name cannot be NULL")

	. <- NULL

	if(family$family=="binomial"){
		if("cbind" %in% all.names(formula)){
		  resp <- paste("cbind(",all.vars(formula)[1],",",all.vars(formula)[2],")")
		  covs <- all.vars(formula)[3:length(all.vars(formula))]
		}
	}else{
		resp <- all.vars(formula)[1]
		covs <- all.vars(formula)[2:length(all.vars(formula))]
	}
  
	bef_df <- sstap_df(stap_formula = stap_formula,
				 subject_data = subject_data,
				 subject_id = subject_id,
				 basis_functions = basis_functions,
				 dt_data = dt_data,
				 BEF_col_name = BEF_col_name,
				 distance_col_name = distance_col_name,
				 time_col_name = time_col_name)

	stap_data <- extract_stap_data(stap_formula)
  
	fixef <- lme4::nobars(formula)
	resp <- all.vars(fixef)[1]
	covs <- all.vars(fixef)[2:length(all.vars(fixef))]
	covs <- c(covs,colnames(bef_df))
	formula <- stats::as.formula(paste(resp, " ~ ", 
							  paste(paste(covs,collapse = " + "), " + ", ## fixef terms
									paste0("(",lme4::findbars(formula),")",collapse=" + ") ))) ## ranef terms

	subject_data <- subject_data %>% dplyr::arrange_(.dots=subject_id)
	X <- cbind(subject_data,bef_df)
	if(method=="lmer"){
		fit <- lme4::lmer(formula,data = X)
	    spatial_covs <- stap_data$covariates[stap_data$stap_code %in% c(0,2)]
	    temp_covs <- stap_data$covariates[stap_data$stap_code %in% c(1,2)]
		fit <- sstapMod(fit, 
					 basis_functions = basis_functions,
					 BEFs = stap_data$covariates,
					 stap_code = stap_data$stap_code,
					 spaceranges = lapply(seq_along(spatial_covs),function(x){dt_data %>% 
					     dplyr::filter(!!dplyr::sym(BEF_col_name) == stap_data$covariates[x]) %>% 
					     dplyr::pull(!!dplyr::sym(distance_col_name)) %>% range(.)}),
					 timeranges = lapply(seq_along(temp_covs),function(x){dt_data %>% 
					     dplyr::filter(!!dplyr::sym(BEF_col_name) == stap_data$covariates[x]) %>% 
					     dplyr::pull(!!dplyr::sym(time_col_name)) %>% range(.)})
					 )
	}
	else if(method=="stan_lmer"){
  
		fit <- rstanarm::stan_lmer(formula,data = X, ...)

		fit$basis_functions <- basis_functions
		fit$stap_data <- extract_stap_data(stap_formula)
		fit$BEFs <- fit$stap_data$covariates
		if(any(fit$stap_data$stap_code %in% c(0,2)))
			fit$spaceranges <- lapply(fit$BEFs,function(x){ dt_data %>% 
			 dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
			 dplyr::pull(!!dplyr::sym(distance_col_name)) %>% range(.)})
		if(any(fit$stap_data$stap_code %in% c(1,2)))
			fit$timeranges <- lapply(fit$BEFs,function(x){ dt_data %>% 
			 dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
			 dplyr::pull(!!dplyr::sym(time_col_name)) %>% range(.)})

		structure(fit, class = c("stanreg","sstap","merMod"))
	}else if(method=="brm"){
		fit <- brms::brm(formula,data = X, ...)

		fit$basis_functions <- basis_functions
		fit$stap_data <- extract_stap_data(stap_formula)
		fit$BEFs <- fit$stap_data$covariates
		if(any(fit$stap_data$stap_code %in% c(0,2)))
			fit$spaceranges <- lapply(fit$BEFs,function(x){ dt_data %>% 
			 dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
			 dplyr::pull(!!dplyr::sym(distance_col_name)) %>% range(.)})
		if(any(fit$stap_data$stap_code %in% c(1,2)))
			fit$timeranges <- lapply(fit$BEFs,function(x){ dt_data %>% 
			 dplyr::filter(!!dplyr::sym(BEF_col_name) == x) %>% 
			 dplyr::pull(!!dplyr::sym(time_col_name)) %>% range(.)})

		structure(fit, class = c("brmsfit","sstap"))
	}
}
