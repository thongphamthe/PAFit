# PAFit package
<!--[![codecov](https://codecov.io/gh/thongphamthe/PAFit/branch/master/graph/badge.svg)](https://codecov.io/gh/thongphamthe/PAFit)-->
[![Downloads from CRAN](https://cranlogs.r-pkg.org/badges/PAFit?color=brightgreen)](https://CRAN.R-project.org/package=PAFit)
[![Downloads from CRAN](https://cranlogs.r-pkg.org/badges/grand-total/PAFit?color=brightgreen)](https://CRAN.R-project.org/package=PAFit)
[![CRAN](https://www.r-pkg.org/badges/version/PAFit)](https://cran.rstudio.com/package=PAFit)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)


This package provides a framework for modelling and inferencing attachment mechanisms of temporal complex networks. For estimating the preferential attachment (PA) function in isolation, we implement Jeong's method, the corrected Newman's method and the PAFit method. For jointly estimating the PA function and node fitnesses, we implement the PAFit method. The package can quantify the remaining uncertainties by providing confidence intervals for the estimated results. We also provide flexible methods to generate a wide range of temporal networks based on PA and fitness.   

Installation
------------

The release version of the package is hosted on [CRAN](https://CRAN.R-project.org/package=PAFit) and can be installed in the usual way:
```r
install.packages("PAFit")
```

This dev version on GitHub can be installed as follows:
```r
require(devtools)
install_github("thongphamthe/PAFit@devel")
```

Getting started
---------------

To get started, load the package
```r
library("PAFit")
```
then work through the tutorial (link to the current CRAN version):

 * [Tutorial](https://CRAN.R-project.org/package=PAFit/vignettes/Tutorial.pdf)
 
An earlier version of this tutorial has been accepted to Journal of Statistical Software:
  * [arXiv version](https://arxiv.org/abs/1704.06017)
  
Reference manual
---------------

Please refer to the current version on CRAN:

 * [Manual](https://CRAN.R-project.org/package=PAFit/PAFit.pdf) 

You can view the html version, which has a better layout but renders mathematical symbols worse than the pdf version, if you use [Rstudio](https://www.rstudio.com/) 

NEWS
---------------

Please refer to the current version on CRAN:

 * [NEWS](https://CRAN.R-project.org/package=PAFit/news.html)

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
