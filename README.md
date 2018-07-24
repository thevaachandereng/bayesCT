# BACT - Bayesian Adaptive Clinical Trial


[![Build Status](https://travis-ci.org/thevaachandereng/BACT.svg?branch=master)](https://travis-ci.org/thevaachandereng/BACT)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://cranlogs.r-pkg.org/badges/BACT)](https://cran.r-project.org/package=BACT)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/thevaachandereng/BACT?branch=master&svg=true)](https://ci.appveyor.com/project/thevaachandereng/BACT)
[![codecov](https://codecov.io/gh/thevaachandereng/BACT/branch/master/graph/badge.svg)](https://codecov.io/gh/thevaachandereng/BACT)
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)


**Authors**: Thevaa Chandereng, Donald Musgrove, Tarek Haddad, Timothy Hanson and Theodore Lystig


Overview
--------

BACT is a R package in Bayesian adaptive clinical trial. 
It is designed to simulate Bayesian adaptive trial for different parameters and distribution.
The BACT website is available [here](https://thevaachandereng.github.io/BACT/). 


Installation
------------
Prior to analyzing your data, the R package needs to be installed.

The easiest way to install BACT is through CRAN:

``` r
install.packages("BACT")
```

There are other additional ways to download BACT.
The first option is most useful if want to download a specific version of BACT
(which can be found at https://github.com/thevaachandereng/BACT/releases).
``` r 
devtools::install_github("thevaachandereng/BACT@vx.xx.x")
# OR 
devtools::install_version("BACT", version = "x.x.x", repos = "http://cran.us.r-project.org")
```

The second option is to download through GitHub. 

``` r
devtools::install_github("thevaachandereng/BACT")
```

After successful installation, the package must be loaded into the working space:

``` r 
library(BACT)
```

Usage
------------
See the vignettes for usage instructions.


License
------------
BACT is available under the open source [GNU General Public License, version 3](https://opensource.org/licenses GPL-3.0).
