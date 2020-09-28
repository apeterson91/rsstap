
## <img src = "docs/figures/rsstap_hex.png" align="right" width="200" height = "200">  `rsstap`: Spline Spatial Temporal Aggregated Predictors
<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Travis build status](https://travis-ci.org/apeterson91/rsstap.svg?branch=master)](https://travis-ci.org/apeterson91/rsstap)
[![R build status](https://github.com/apeterson91/rsstap/workflows/R-CMD-check/badge.svg)](https://github.com/apeterson91/rsstap/actions)
<!-- badges: end -->

## About

This is an R package that fits spline spatial temporal aggregated predictors to subject indexed spatio-temporal data.
The primary target audience is researchers interested in the effect of built environment features (BEFs) on human health,though other applications are possible.

## Installation

### Development Version

As this package is under active development it is currently only available via Github. If you'd like to to install the current development software use the following 
 lines of R code

 ```r
 if(!require(devtools)){
	install.packages("devtools")
	library(devtools)
 }

install_github("apeterson91/rsstap",dependencies = TRUE)
 ```


#### Code of Conduct

Please note that `rsstap` is released with a [Contributor Code of Conduct](https://www.contributor-covenant.org/). By contributing to this project, you agree to abide by its terms.


## How to cite this package

 A citation is in progress. Check back soon.

## Acknowledgments

This work was developed with support from NIH grant R01-HL131610 (PI: Sanchez).


