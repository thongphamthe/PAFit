\name{PAFit-package}
\alias{PAFit-package}
\alias{summary.CV_Data}
\alias{summary.PAFit_data}
\alias{summary.PA_result}
\alias{summary.PAFit_result}
\alias{summary.CV_Result}
\alias{print.CV_Data}
\alias{print.PAFit_data}
\alias{print.PA_result}
\alias{print.PAFit_result}
\alias{print.CV_Result}
\alias{PAFit}
\alias{PAFit_new}
\docType{package} 
\title{Joint Inference of Preferential Attachment and Node Fitness in Temporal Complex Networks}
\description{
An implementation of a framework for modelling and inferencing the attachment mechanisms of temporal complex networks. For estimating the preferential attachment (PA) function in isolation, we implement Jeong's method, the corrected Newman's method and the PAFit method. For jointly estimating the PA function and node fitnesses, we implement the PAFit method. The package also provides flexible methods to generate a wide range of temporal networks based on PA and fitness. See the accompanying vignette for a tutorial. For a list of references, please run the command: \code{citation("PAFit")}.
}
\details{
\tabular{ll}{
Package: \tab PAFit\cr
Type: \tab Package\cr
Version: \tab 0.9.7.5\cr
Date: \tab 2017-03-01\cr
License: \tab GPL-3 \cr
}
The main functions are as follows. 

Basic utilities:
\itemize{
\item \code{\link{GenerateNet}}: generates a wide range of networks based on the preferential attachment and fitness mechanisms.
\item \code{\link{GetStatistics}}: summarizes a matrix of edges into summary statistics ready for estimating the preferential attachment function or node fitness
}
Functions to estimate the attachment function in isolation:
\itemize{
\item \code{\link{Jeong}}: estimates the preferential attachment in isolation by Jeong's method (Ref. 1)
\item \code{\link{Newman}}: estimates the preferential attachment in isolation by the corrected Newman's method (Refs. 2, 3)
\item \code{\link{OnlyA_Estimate}}: estimates the preferential attachment in isolation by PAFit method (Ref. 3)
}
Function to estimate node fitnesses in isolation:
\itemize{
\item \code{\link{OnlyF_Estimate}}: estimates the preferential attachment in isolation by Jeong's method (Ref. 4)
}
Function to estimate the attachment function and node fitnesses jointly:
\itemize{
\item \code{\link{JointEstimate}}: jointly estimates the preferential attachment function and node fitness by PAFit method (Ref. 4)
}

}

\author{
Thong Pham, Paul Sheridan, Hidetoshi Shimodaira.
Maintainer: Thong Pham \email{thongpham@thongpham.net}
}


\references{
1. Jeong, H., \enc{Néda}{Neda}, Z. & \enc{Barabási}{Barabasi}, A. . Measuring preferential attachment in evolving networks. Europhysics Letters. 2003;61(61):567–572. doi: 10.1209/epl/i2003-00166-9 (\url{http://iopscience.iop.org/article/10.1209/epl/i2003-00166-9/fulltext/}).

2. Newman, M.. Clustering and preferential attachment in growing networks. Physical Review E. 2001;64(2):025102 (\url{https://journals.aps.org/pre/abstract/10.1103/PhysRevE.64.025102}).

3. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).

4. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{www.nature.com/articles/srep32558}).
}


\concept{temporal complex networks}
\concept{preferential attachment}
\concept{attachment function}
\concept{fitness model}
\concept{Jeong's method}
\concept{Corrected Newman's method}
\concept{\enc{Barabási}{Barabasi}-Albert model}
\concept{Bianconi-\enc{Barabási}{Barabasi} model}

\examples{
\dontrun{
### Jointly estimate the attachment function and node fitnesses
library("PAFit")
  net        <- GenerateNet(N = 1000 , m = 10 , mode = 1 , alpha = 1 , shape = 5, rate = 5)
  net_stats  <- GetStatistics(net$graph)
  result     <- JointEstimate(net$graph, net_stats)
  
   # true function
  true_A     <- pmax(result$estimate_result$center_k,1)
  #plot the estimated PA function
  plot(result$estimate_result , net_stats)
  lines(result$estimate_result$center_k + 1, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
  
  #plot the estimated node fitnesses and true node fitnesses
  plot(result$estimate_result , net_stats, true = net$fitness, plot = "true_f")
}
}
