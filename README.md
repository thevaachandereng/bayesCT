bayesCT - Tool for Simulation and Analysis of Adaptive Bayesian Clinical Trials <img src="man/figures/logo.png" align="right" width="150" height="150" />
===============================================================================


[![Build Status](https://travis-ci.org/thevaachandereng/bayesCT.svg?branch=master)](https://travis-ci.org/thevaachandereng/bayesCT)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/bayesCT)](https://cran.r-project.org/package=bayesCT)
[![Download_Badge](https://cranlogs.r-pkg.org/badges/bayesCT)](https://cran.r-project.org/package=bayesCT)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/donaldmusgrove/bayesDP/issues)
[![Build status](https://ci.appveyor.com/api/projects/status/2wfwigrrcpom0oi9/branch/master?svg=true)](https://ci.appveyor.com/project/thevaachandereng/bayesct/branch/master)
[![codecov](https://codecov.io/gh/thevaachandereng/bayesCT/branch/master/graph/badge.svg)](https://codecov.io/gh/thevaachandereng/bayesCT)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/thevaachandereng/bayesCT/master?urlpath=rstudio)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)


**Authors**: Thevaa Chandereng, Donald Musgrove, Tarek Haddad, Graeme Hickey, Timothy Hanson and Theodore Lystig


Overview
--------
`bayesCT` is a R package for simulation and analysis of adaptive Bayesian randomized controlled trials under a range of trial designs and outcome types. Currently, it supports Gaussian, binomial, and time-to-event. The `bayesCT` package website is available [here](https://thevaachandereng.github.io/bayesCT/). Please note this package is still under development. 


Installation
------------
Prior to analyzing your data, the R package needs to be installed. The easiest way to install `bayesCT` is through CRAN:

``` r
install.packages("bayesCT")
```

There are other additional ways to download `bayesCT`. The first option is most useful if want to download a specific version of `bayesCT` (which can be found at https://github.com/thevaachandereng/bayesCT/releases):

``` r 
devtools::install_github("thevaachandereng/bayesCT@vx.xx.x")

# or 

devtools::install_version("bayesCT", version = "x.x.x", repos = "http://cran.us.r-project.org")
```

The second option is to download through GitHub: 

``` r
devtools::install_github("thevaachandereng/bayesCT")
```

After successful installation, the package must be loaded into the working space:

``` r 
library(bayesCT)
```


Usage
------------
See the vignettes for usage instructions and example.


Code of Conduct
------------
Please note that this project is released with a [Contributor Code of Conduct](https://github.com/thevaachandereng/bayesCT/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.


License
------------
`bayesCT` is available under the open source [GNU General Public License, version 3](https://www.r-project.org/Licenses/GPL-3).
