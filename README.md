# The PAFit package
[![Build Status](https://travis-ci.org/thongphamthe/PAFit.png?branch=master)](https://travis-ci.org/thongphamthe/PAFit)
[![codecov](https://codecov.io/gh/thongphamthe/PAFit/branch/master/graph/badge.svg)](https://codecov.io/gh/thongphamthe/PAFit)
[![Downloads](http://cranlogs.r-pkg.org/badges/PAFit?color=brightgreen)](http://cran.rstudio.com/package=PAFit)
[![CRAN](http://www.r-pkg.org/badges/version/PAFit)](http://cran.rstudio.com/package=PAFit)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

====================

This package provides a framework for modelling and inferencing the attachment mechanisms of temporal complex networks. For estimating the preferential attachment (PA) function in isolation, it implements Jeong's method, the corrected Newman's method and the PAFit method. For jointly estimating the PA function and node fitnesses, it implements the PAFit method. It also provides flexible methods to generate a wide range of temporal networks based on PA and fitness.   

Installation
------------

This package is hosted on [CRAN](http://cran.r-project.org/web/packages/PAFit/) and can be installed in the usual way:
```r
install.packages("PAFit")
```
Alternatively, the development version can be install from from github using the devtools package:
```r
install.packages("devtools")
devtools::install_github("thongphamthe/PAFit", sub = "pkg")
```

Note: Windows users have to first install [Rtools](http://cran.rstudio.com/bin/windows/Rtools/).

Getting Started
---------------

To get started, load the package
```r
library("PAFit")
```
then work through the tutorial (link to the current CRAN version):

 * [Tutorial](https://cran.r-project.org/web/packages/PAFit/vignettes/Tutorial.pdf)

NEWS
---------------

Please refer to the current CRAN version for news:

 * [NEWS](https://cran.r-project.org/web/packages/PAFit/NEWS)

Citation
---------------

Please use the information from the citation information file (link to the current CRAN version): 

 * [Citation Info](https://cran.r-project.org/web/packages/PAFit/citation.html)
 
License
-----------------
GPL-3

Other information
-----------------

 * If you have any suggestions or find bugs, please use the github [issue tracker](https://github.com/thongphamthe/PAFit/issues)
 * Feel free to submit pull requests
