# Part of the rsstap package for estimating model parameters
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
#
#' Datasets for rstap examples
#' 
#' Small datasets for use in \pkg{rstap} examples and vignettes.
#'
#' @name rsstap-datasets 
#' @aliases rsstap_demo rsstap_demo rsstap_repeatedmeasures 
#'
#' @format 
#' \describe{
#' \item{\code{rsstap_demo}}{
#'	\itemize{
#'          \item \code{subject_data}
#'			\item -  \code{subj_id}: The subject unique identifier
#'			\item -  \code{sex}: Subject binary covariate
#'			\item -  \code{y}: Subject outcome
#'	        \item \code{direct_distance_data}
#'			\item -  \code{subj_id}: Subject unique identifier
#'			\item -  \code{BEF}: BEF categorical coding
#'			\item -  \code{Distance}: Distance between subject and BEF
#'		\item Miscl
#'		\item -  alpha: simulated intercept
#'		\item -  delta: simulated binary covariate effect
#'		\item -  sigma: simulated standard deviation for residual normal errors
#'		\item -  true effect: simulated true effect
#'		}
#'	}
#' Source: \href{https://biostatistics4socialimpact.github.io/rsstap/articles/Introduction}{Introductory Vignette}
#' 
#' \item{\code{rsstap_repeatedmeasures}}{
#' \itemize{
#'	\item Still in Progress
#' }}
#'}
#' 
#'
NULL
