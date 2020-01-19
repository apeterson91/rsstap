# bbnet: Bayesian Built Environment Network Models
 <!--badges: start -->
<!-- Travis badge: start  [![Travis build status](https://travis-ci.org/apeterson91/bendr.svg?branch=master)](https://travis-ci.org/apeterson91/bendr)-->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## About

This is an R package that fits the bayesian built environment network model to
subject indexed spatio-temporal data. The primary target audience is researchers interested in the effect of built environment features (BEFs) on human health, 
though other applications are possible.

## Installation

### Development Version

 Currently this package is only available via Github and is not ready for regular use. If you'd like to to install the current development software use the following 
 lines of R code

 ```r
 if(!require(devtools)){
	install.packages("devtools")
	library(devtools)
 }

install_github("apeterson91/bbnet",dependencies = TRUE)
 ```


#### Code of Conduct

Please note that `bbnet` is released with a [Contributor Code of Conduct](https://www.contributor-covenant.org/). By contributing to this project, you agree to abide by its terms.


## How to cite this package

 A citation is in progress. Check back soon.

## Acknowledgments

This work was developed with support from NIH grant R01-HL131610 (PI: Sanchez).


