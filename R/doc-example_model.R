# Part of the rsstap package for investigating features of 
# the built environment
# Copyright (C) 2020 Adam Peterson

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Example model
#' 
#' A model for use in \pkg{rsstap} examples. 
#' 
#' @name example_model
#' @format Calling \code{example("example_model")} will run the model in the 
#'   Examples section, below, and the resulting sstapreg  object will then be
#'   available in the global environment. The \code{chains} and \code{iter}
#'   arguments are specified to make this example be small in size. In practice,
#'   it is reccomended that they be left unspecified in order to use the default
#'   values (4 and 2000 respectively) or increased if there are convergence
#'   problems. The \code{cores} argument is optional and on a multicore system,
#'   the user may well want to set that equal to the number of chains being
#'   executed.
#'   
#' @seealso \code{\link[rbenvo]{example_benvo}} for a description of the data.
#' @examples
#' example_model <- 
#'   sstap_lm(BMI ~ sex + sap(FFR),
#'              benvo = rbenvo::FFbenvo,
#'              # this next line is only to keep the example small in size!
#'              chains = 1, cores = 1, seed = 12345, iter = 500, refresh = 0)
#' example_model
NULL
