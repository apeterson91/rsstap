# Part of the bbnet package 
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#' Create a stap object 
#'
#' @param object  a named list of objects containing information about the staps for a given model 
#' @return an object of class "stap"
stap <- function(object) {
    
    stap_code <- array(sapply(object, function(x) x$stap_code), dim = length(object))
	Q <- length(stap_code)
    covariates <- sapply(object, function(x) x$covariate) 
    dnd_code <- array(sapply(object, function(x) x$dnd_code), dim = Q)
    bar_code <- array(sapply(object, function(x) x$bar_code), dim = Q)
	group_indicator <- sapply(object, function(x) x$group_var)
	group_term <- sapply(object,function(x) x$group_term)
    out <- list(covariates = covariates,
                 stap_code = stap_code,
                 dnd_code = dnd_code,
                 bar_code = bar_code,
                group_indicator = group_indicator,
                group_terms = group_term
	)
    out <- structure(out, class = c("stap"))
    check_dups(out)
    return(out)
}

sap_covs <- function(x)
    UseMethod("sap_covs")

tap_covs <- function(x)
    UseMethod("tap_covs")

stap_covs <- function(x)
    UseMethod("stap_covs")

check_dups <- function(x)
	UseMethod("check_dups")

sap_covs.stap <- function(x)
	return(x$covariates[intersect(which(x$stap_code==0),which(x$group_indicator==0))])

stap_covs.stap <- function(x)
	return(x$covariates[intersect(which(x$stap_code==2),which(x$group_indicator==0))])

tap_covs.stap <- function(x)
	return(x$covariates[intersect(which(x$stap_code==1),which(x$group_indicator==0))])


check_dups.stap <- function(object){
    sap <- sap_covs(object)
    tap <- tap_covs(object)
    stap <- stap_covs(object)
    if(length(sap) != length(unique(sap)))
        stop("A BEF may only be given one kind of stap specification
             either (exclusively) stap or tap or sap")
    else if(length(tap) != length(unique(tap)))
        stop("A BEF may only be given one kind of stap specification
             either (exclusively) stap or tap or sap")
    else if(length(stap) != length(unique(stap)))
       stop("A BEF may only be given one kind of stap specification
             either (exclusively) stap or tap or sap")
    else if(!all.equal(union(sap,stap),c(sap,stap)))
       stop("A BEF may only be given one kind of stap specification
             either (exclusively) stap or tap or sap")
    else if(!all.equal(union(tap,stap),c(tap,stap)))
        stop("A BEF may only be given one kind of stap specification
             either (exclusively) stap or tap or sap")
}
