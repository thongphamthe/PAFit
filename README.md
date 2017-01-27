# The PAFit package
[![Build Status](https://travis-ci.org/thongphamthe/PAFit.png?branch=master)](https://travis-ci.org/thongphamthe/PAFit)
[![codecov](https://codecov.io/gh/thongphamthe/PAFit/branch/master/graph/badge.svg)](https://codecov.io/gh/thongphamthe/PAFit)
[![Downloads](http://cranlogs.r-pkg.org/badges/PAFit?color=brightgreen)](http://cran.rstudio.com/package=PAFit)
[![CRAN](http://www.r-pkg.org/badges/version/PAFit)](http://cran.rstudio.com/package=PAFit)

====================

This package implements popular methods to estimate the preferential attachment function and node fitnesses in a growing network.
In particular, it implements Jeong's method, Newman's method and PAFit.

Installation
------------

This package is hosted on [CRAN](http://cran.r-project.org/web/packages/PAFit/) and can be installed in the usual way:
```r
install.packages("PAFit")
```
Alternatively, the development version can be install from from github using the devtools package:
```r
install.packages("devtools")
devtools::install_github("thongphamthe/PAFit", sub="pkg")
```

Note Windows users have to first install [Rtools](http://cran.rstudio.com/bin/windows/Rtools/).

Getting Started
---------------

To get started, load the package
```r
library("PAFit")
```
then work through the tutorial (links to the current CRAN version):

 * [Tutorial](https://cran.r-project.org/web/packages/PAFit/vignettes/Tutorial.pdf)

License
-----------------
GPL-3

Other information
-----------------

 * If you have any suggestions or find bugs, please use the github [issue tracker](https://github.com/thongphamthe/PAFit/issues)
 * Feel free to submit pull requests
