% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_xsection}
\alias{plot_xsection}
\alias{plot_xsection.sstapreg}
\title{Plot Cross-Sections}
\usage{
plot_xsection(
  x,
  stap_term = NULL,
  component = "Distance",
  fixed_val = 1,
  p = 0.95,
  grid = NULL
)

\method{plot_xsection}{sstapreg}(
  x,
  stap_term = NULL,
  component = "Distance",
  fixed_val = 1,
  p = 0.95,
  grid = NULL
)
}
\arguments{
\item{x}{a sstapreg object}

\item{stap_term}{name of stap term to plot}

\item{component}{one of c("Distance","Time")}

\item{fixed_val}{vector that contains fixed values for whichever component was not specified}

\item{p}{probability_interval}

\item{optional}{grid by default is NULL corresponding to a .1 step grid from the min-max of the covariate space}
}
\description{
Plot Cross-Sections

Plot Cross-Sections
}
\section{Methods (by class)}{
\itemize{
\item \code{sstapreg}: 
}}

\keyword{internal}
