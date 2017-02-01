# PAFit package
[![Build Status](https://travis-ci.org/thongphamthe/PAFit.svg?branch=devel)](https://travis-ci.org/thongphamthe/PAFit)
[![Build status](https://ci.appveyor.com/api/projects/status/ufje8pyddp42tbfu/branch/master?svg=true)](https://ci.appveyor.com/project/thongphamthe/pafit/branch/devel)
[![codecov](https://codecov.io/gh/thongphamthe/PAFit/branch/master/graph/badge.svg)](https://codecov.io/gh/thongphamthe/PAFit)
[![Downloads from CRAN](https://cranlogs.r-pkg.org/badges/PAFit?color=brightgreen)](https://CRAN.R-project.org/package=PAFit)
[![CRAN](https://www.r-pkg.org/badges/version/PAFit)](https://cran.rstudio.com/package=PAFit)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)

====================

This package provides a framework for modelling and inferencing attachment mechanisms of temporal complex networks. For estimating the preferential attachment (PA) function in isolation, we implement Jeong's method, the corrected Newman's method and the PAFit method. For jointly estimating the PA function and node fitnesses, we implement the PAFit method. The package can quantify the remaining uncertainties by providing confidence intervals for the estimated results. We also provide flexible methods to generate a wide range of temporal networks based on PA and fitness.   

Installation
------------

This package is hosted on [CRAN](https://CRAN.R-project.org/package=PAFit) and can be installed in the usual way:
```r
install.packages("PAFit")
```
Alternatively, the development version can be installed from github using the devtools package (Windows users have to first install [Rtools](https://cran.rstudio.com/bin/windows/Rtools/)):
```r
install.packages("devtools")
devtools::install_github("thongphamthe/PAFit")
```

Getting started
---------------

To get started, load the package
```r
library("PAFit")
```
then work through the tutorial (link to the current CRAN version):

 * [Tutorial](https://CRAN.R-project.org/package=PAFit/vignettes/Tutorial.pdf)
 
Reference manual
---------------

Please refer to the current version on CRAN:

 * [Manual](https://CRAN.R-project.org/package=PAFit/PAFit.pdf) 

NEWS
---------------

Please refer to the current version on CRAN:

 * [NEWS](https://CRAN.R-project.org/package=PAFit/NEWS)

Citation
---------------

Please refer to the citation information file (link to the current CRAN version): 

 * [Citation Info](https://CRAN.R-project.org/package=PAFit/citation.html)
 
License
-----------------
GPL-3

Other information
-----------------

 * If you have any suggestions or find bugs, please use the github [issue tracker](https://github.com/thongphamthe/PAFit/issues)
 * Feel free to submit pull requests
