\name{PAFit-package}
\alias{PAFit-package}
\docType{package} 
\title{Joint Inference of Preferential Attachment and Node Fitness in Temporal Complex Networks}
\description{
An implementation of a framework for modelling and inferencing the attachment mechanisms of temporal complex networks. For estimating the preferential attachment (PA) function in isolation, we implement Jeong's method, the corrected Newman's method and the PAFit method. For jointly estimating the PA function and node fitnesses, we implement the PAFit method. The package also provides flexible methods to generate a wide range of temporal networks based on PA and fitness. See the accompanying vignette for a tutorial. For a list of references, please run the command: \code{citation("PAFit")}.
}
\details{
\tabular{ll}{
Package: \tab PAFit\cr
Type: \tab Package\cr
Version: \tab 0.9.7\cr
Date: \tab 2017-02-15\cr
License: \tab GPL-3 \cr
}
The main functions are:
\itemize{
\item \code{\link{GenerateNet}}: generates a wide range of networks based on the preferential attachment and fitness mechanisms.
\item \code{\link{GetStatistics}}: summarizes a matrix of edges into summary statistics ready for estimating the preferential attachment function or node fitness
\item \code{\link{Jeong}}: estimates the preferential attachment in isolation by Jeong's method 
\item \code{\link{Newman_corrected}}: estimates the preferential attachment in isolation by the corrected Newman's method
\item \code{\link{PAFit}}: jointly estimates the preferential attachment function and node fitness
\item \code{\link{CreateDataCV}}: creates cross-validation for \code{\link{performCV}} function
\item \code{\link{performCV}}: performs cross-validation for selecting regularization parameters of \code{\link{PAFit}} function
}
}
\author{
Thong Pham, Paul Sheridan, Hidetoshi Shimodaira.
Maintainer: Thong Pham \email{thongpham@thongpham.net}
}


\references{
1. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Nonparametric Estimation of the Preferential Attachment Function in Complex Networks: Evidence of Deviations from Log Linearity, Proceedings of ECCS 2014, 141-153 (Springer International Publishing) (\url{http://dx.doi.org/10.1007/978-3-319-29228-1_13}).

2. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).

3. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{www.nature.com/articles/srep32558}).
}


\keyword{temporal complex networks}
\keyword{preferential attachment}
\keyword{fitness model}
\keyword{Jeong's method}
\keyword{Corrected Newman's method}
\keyword{\enc{Barabási}{Barabasi}-Albert model}
\keyword{Bianconi-\enc{Barabási}{Barabasi} model}

\examples{
library("PAFit")
net        <- GenerateNet(N = 50,m = 5,alpha = 0.5, mode = 1, shape = 100, rate = 100)
net_stats  <- GetStatistics(net$graph)
result     <- PAFit(net_stats, s = 100)
summary(result)
plot(result,net_stats)
}
