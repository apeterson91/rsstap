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

#' @title Lmer rsstap model objects
#'
#' @description
#' A fitted lmer rsstap object. Inherits from \pkg{lme4}'s
#' \code{\link[lme4]{merMod}} so that most methods defined for this
#' parent class should work seamlessly within \pkg{rsstap} analysis
#'
#' @slot basis_functions 
#'      list of basis functions used for each STAP
#' @slot BEFs
#'   A character vector containing the names of the BEFs included in the model fit
#' @slot stap_code
#'   stap_code designation
#' @slot spaceranges 
#'    list of spatial range for each STAP
#' @slot timeranges 
#'    list of temporal ranges for each STAP 
#'
#' @export
#'
sstapMod <- setClass("sstapMod",
                  slots = c(
                    basis_functions = "list",
                    BEFs = "character",
				    stap_code = "array",
					spaceranges = "list",
					timeranges = "list"
                  ),
                  contains = "lmerMod"  # also from lme4
				)
