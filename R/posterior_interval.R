# This software is part of the rsstap package
# Copyright (C) 2020 Adam Peterson
#
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

#' Posterior uncertainty intervals
#'
#' The \code{posterior_interval} function computes Bayesian posterior uncertainty
#' intervals. These intervals are often referred to as \emph{credible}
#' intervals. This text is inspired by the same function documentation from \pkg{rstanarm}.
#'
#' @export
#'
#' @param object sstapreg object
#' @param prob A number \eqn{p \in (0,1)}{p (0 < p < 1)} indicating the desired
#'   probability mass to include in the intervals. The default is to report
#'   \eqn{95}\% intervals (\code{prob=0.95}).
#' @param pars vector of parameter names
#' @param ... ignored
#'
#' @return A matrix with two columns and as many rows as model parameters (or
#'   the subset of parameters specified by \code{pars}. 
#'   For a given value of \code{prob}, \eqn{p}, the columns
#'   correspond to the lower and upper \eqn{100p}\% interval limits and have the
#'   names \eqn{100\alpha/2}\% and \eqn{100(1 - \alpha/2)}\%, where \eqn{\alpha
#'   = 1-p}. For example, if \code{prob=0.9} is specified (a \eqn{90}\%
#'   interval), then the column names will be \code{"5\%"} and \code{"95\%"},
#'   respectively.
#'
#' @details
#' \subsection{Interpretation}{
#' Unlike for a frenquentist confidence interval, it is valid to say that,
#' conditional on the data and model, we believe that with probability \eqn{p}
#' the value of a parameter is in its \eqn{100p}\% posterior interval. This
#' intuitive interpretation of Bayesian intervals is often erroneously applied
#' to frequentist confidence intervals. See Morey et al. (2015) for more details
#' on this issue and the advantages of using Bayesian posterior uncertainty
#' intervals (also known as credible intervals).
#' }
#'
#'
#' @template reference-gelman-carlin
#' @template reference-morey
#' @importFrom rstantools posterior_interval
#'
posterior_interval.sstapreg <-
  function(object,
           prob = 0.95,
           pars = NULL,
           ...) {
    
    mat <- as.matrix.sstapreg(object, pars = pars)
    rstantools::posterior_interval(mat, prob = prob)
  }

#' @export
rstantools::posterior_interval
