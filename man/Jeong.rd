\name{Jeong}
\alias{Jeong}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Jeong's method for estimating  the preferential attachment function
}
\description{
  This function estimates the preferential attachment function by Jeong's method. 
}
\usage{
    Jeong(raw_net, net_stat, T_0, T_1 , interpolate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw_net}{
    a three-column matrix that contains the network.
  }
  \item{net_stat}{
    An object of class "PAFit_data" which contains summerized statistics needed in estimation.
  }
  \item{T_0}{Positive integer. The starting time-step at which the degree sequence is freezed.}
  \item{T_1}{Positive integer. The ending time-step.
  } 
  \item{interpolate}{
    Logical. If \code{TRUE} then all the gaps in the estimated PA function are interpolated by linear interpolating in logarithm scale. Default value is \code{FALSE}.
  }
}
\value{
  Outputs an "PA_result" object which contains the estimated PA function.
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Nonparametric Estimation of the Preferential Attachment Function in Complex Networks: Evidence of Deviations from Log Linearity, Proceedings of ECCS 2014, 141-153 (Springer International Publishing) (\url{http://dx.doi.org/10.1007/978-3-319-29228-1_13}).
  
  2. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).
  
  3. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{www.nature.com/articles/srep32558}).
}

\examples{
  library("PAFit")
  net        <- GenerateNet(N = 1000 , m = 1 , mode = 1 , alpha = 1 , shape = 0)
  net_stats  <- GetStatistics(net$graph)
  result     <- Jeong(net$graph, net_stats, T_0 = 500, T_1 = 900)
  summary(result)
}
