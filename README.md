# bayesCT - Tool for Simulation in Adaptive Bayesian Clinical Trials


[![Build Status](https://travis-ci.org/thevaachandereng/bayesCT.svg?branch=master)](https://travis-ci.org/thevaachandereng/bayesCT)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://cranlogs.r-pkg.org/badges/bayesCT)](https://cran.r-project.org/package=bayesCT)
[![Build status](https://ci.appveyor.com/api/projects/status/2wfwigrrcpom0oi9/branch/master?svg=true)](https://ci.appveyor.com/project/thevaachandereng/bayesct/branch/master)
[![codecov](https://codecov.io/gh/thevaachandereng/bayesCT/branch/master/graph/badge.svg)](https://codecov.io/gh/thevaachandereng/bayesCT)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)


**Authors**: Thevaa Chandereng, Donald Musgrove, Tarek Haddad, Timothy Hanson and Theodore Lystig


Overview
--------

bayesCT is a R package for simulation in Adaptive Bayesian Clinical Trials
It is designed to simulate Bayesian adaptive trial for different parameters and distribution.
Currently, it supports normal, binomial and survival data.
The BACT website is available [here](https://thevaachandereng.github.io/bayesCT/). 
The package is still under development. 


Installation
------------
Prior to analyzing your data, the R package needs to be installed.

The easiest way to install bayesCT is through CRAN:

``` r
install.packages("bayesCT")
```

There are other additional ways to download bayesCT.
The first option is most useful if want to download a specific version of BACT
(which can be found at https://github.com/thevaachandereng/bayesCT/releases).
``` r 
devtools::install_github("thevaachandereng/bayesCT@vx.xx.x")
# OR 
devtools::install_version("bayesCT", version = "x.x.x", repos = "http://cran.us.r-project.org")
```

The second option is to download through GitHub. 

``` r
devtools::install_github("thevaachandereng/bayesCT")
```

After successful installation, the package must be loaded into the working space:

``` r 
library(bayesCT)
```

Usage
------------
See the vignettes for usage instructions.


License
------------
bayesCT is available under the open source [GNU General Public License, version 3](https://www.r-project.org/Licenses/GPL-3).
