\name{Newman_corrected}
\alias{Newman_corrected}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Corrected Newman's method for estimating the preferential attachment function
}
\description{
This function estimates the preferential attachment function by the corrected Newman's method. 
}
\usage{
  Newman_corrected(net_stat , start = 1 , interpolate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net_stat}{
    An object of class \code{PAFit_data} which contains summerized statistics needed in estimation. This object is created by the function \code{\link{GetStatistics}}.
  }
  \item{start}{Positive integer. The starting time from which the method is applied. Default value is \eqn{1}.}
  \item{interpolate}{
    Logical. If \code{TRUE} then all the gaps in the estimated PA function are interpolated by linear interpolating in logarithm scale. Default value is \code{FALSE}.
  }
}
\value{
  Outputs an \code{PA_result} object which contains the estimated PA function. It also includes the estimated attachment exponenent \eqn{\alpha} and the confidence interval of \eqn{\alpha} when possible.
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Nonparametric Estimation of the Preferential Attachment Function in Complex Networks: Evidence of Deviations from Log Linearity, Proceedings of ECCS 2014, 141-153 (Springer International Publishing) (\url{http://dx.doi.org/10.1007/978-3-319-29228-1_13}).
  
  2. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).
  
  3. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{www.nature.com/articles/srep32558}).
}
\seealso{\code{\link{Jeong}}, \code{\link{PAFit}}}
\examples{
  library("PAFit")
  net        <- GenerateNet(N = 1000 , m = 30 , mode = 1 , alpha = 1 , shape = 0) # no fitness
  net_stats  <- GetStatistics(net$graph)
  result     <- Newman_corrected(net_stats)
  summary(result)
}
