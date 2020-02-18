# extract stap data from formula and create stap object 
# 
# @param formula that designates model expression including stap covariates 
#
extract_stap_data <- function(formula){

	  nbf <- lme4::nobars(formula)
    all_names <- all.names(nbf)
    staps <- c("sap","tap","stap")
    staps <- c(staps,paste0(staps,"_log"),paste0(staps,"_dnd"),paste0(staps,"_bar"),
               paste0(staps,"_dnd_bar"),paste0(staps,"_bar_dnd"))
    stap_covs <- all_names[which(all_names %in% staps) +1]
    stap_code <- get_stap_code(all_names,stap_covs)
    dnd_code <- sapply(all_names[which(all_names %in% staps)], function(x) grepl("_dnd",x)*1)
    bar_code <- sapply(all_names[which(all_names %in% staps)], function(x) grepl("_bar",x)*1)
  	wbf <- lme4::findbars(formula)
  	all_names <- unlist(lapply(wbf,all.names))
  	stap_group_covs <- all_names[which(all_names %in% staps) +1]
    if(length(stap_covs)==0 & length(stap_group_covs)==0)
        stop("No stap covariates specified")
	stap_group_term <- all_names[which(all_names %in% staps) + 2]
	stap_group_code <- get_stap_code(all_names,stap_group_covs)
    out <- lapply(seq_along(stap_covs),function(x){ list(covariate = stap_covs[x],
                                                    stap_code = stap_code[x],
                                                    dnd_code = dnd_code[x],
                                                    bar_code = bar_code[x],
													group_var = 0,
													group_term = NA)})
    if(length(out)>0)
  	  out <- c(out,lapply(seq_along(stap_group_covs),function(x){ list(covariate = stap_group_covs[x],
  																	stap_code = stap_group_code[x],
  																	group_var = 1,
  																	dnd_code = 0,
  																	bar_code = 0,
  																	group_term = stap_group_term[x]
  																	)}))
    else
      out <- lapply(seq_along(stap_group_covs),function(x){ list(covariate = stap_group_covs[x],
                                                                       stap_code = stap_group_code[x],
                                                                       group_var = 1,
                                                                       dnd_code = 0,
                                                                       bar_code = 0,
                                                                       group_term = stap_group_term[x])})
      
	  out <- stap(out)
	  check_dups(out)
    return(out)
}

# Get stap coding from formula
#
# @param  all_names character vector from calling all.names(formula)
# @param names of the stap covariates
# @return vector of length equal to number of staps + saps + taps
# with the appropriate coding for each appropriate predictor
get_stap_code <- function(all_names,stap_covs){

    staps <- c("sap" = 0,"tap" = 1,"stap" = 2,
               "sap_log" = 0, "tap_log" = 1,
				"stap_log" = 2,"sap_dnd" = 0,
				"tap_dnd" = 1, "stap_dnd" = 2,
                "sap_bar" = 0, "tap_bar" = 1,
				"stap_bar" = 2,"sap_dnd_bar" = 0,
				"tap_dnd_bar" = 1, "stap_dnd_bar" = 2,
               "sap_bar_dnd" = 0, "tap_bar_dnd" = 1,
			   "stap_dnd_bar" = 2)
    sapply(unique(stap_covs),function(x) as.vector(staps[all_names[which(all_names == x)-1]]))
}


create_stap_lmer_formula <- function(stap_formula,BEF_names,stap_group_var,stap_group_terms,coef_names){
  if(length(lme4::findbars(stap_formula))==0)
    return(list())
  
	stap_group_ics <- which(stap_group_var==1)
	stap_group_covs <- paste0(stap_group_terms[stap_group_ics],"_",BEF_names[stap_group_ics],"_[0-9]")
	BEF_cov_ics <- lapply(stap_group_covs,function(x) grep(x,coef_names))
	formula_term <- lapply(1:length(BEF_cov_ics),function(x) paste0("(",paste0(coef_names[BEF_cov_ics[[x]]], collapse = " + ") ,"|",stap_group_terms[stap_group_ics][x], " )"))
	formula_term <- stringr::str_c(unlist(formula_term),collapse=" + ")
	return(formula_term)
}
