\name{Newman}
\alias{Newman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Corrected Newman's method for estimating the preferential attachment function
}
\description{
This function implements a correction of Newman's method to estimate the preferential attachment function. 
}
\usage{
  Newman(raw_net, net_stat , start = 1 , interpolate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw_net}{
    a three-column matrix that contains the network.
}
  \item{net_stat}{
    An object of class \code{PAFit_data} which contains summerized statistics needed in estimation. This object is created by the function \code{\link{GetStatistics}}.
  }
  \item{start}{Positive integer. The starting time from which the method is applied. Default value is \eqn{1}.}
  \item{interpolate}{
    Logical. If \code{TRUE} then all the gaps in the estimated PA function are interpolated by linear interpolating in logarithm scale. Default value is \code{FALSE}.
  }
}
\value{
  Outputs an \code{PA_result} object which contains the estimated attachment function. It also includes the estimated attachment exponenent \eqn{\alpha} (the field \code{alpha}) and the confidence interval of \eqn{\alpha} (the field \code{ci}) when possible.
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Newman, M.. Clustering and preferential attachment in growing networks. Physical Review E. 2001;64(2):025102 (\url{https://journals.aps.org/pre/abstract/10.1103/PhysRevE.64.025102}).
}
\seealso{

See \code{\link{GetStatistics}} for how to create summerized statistics needed in this function.
 
See \code{\link{Jeong}}, \code{\link{OnlyA_Estimate}} for other methods to estimate the attachment function in isolation.
}
\examples{
  library("PAFit")
  net        <- GenerateNet(N = 1000 , m = 1 , mode = 1 , alpha = 1 , shape = 0)
  net_stats  <- GetStatistics(net$graph)
  result     <- Newman(net$graph,net_stats)
  # true function
  true_A     <- result$center_k
  #plot the estimated PA function
  plot(result , net_stats)
  lines(result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
}

\concept{preferential attachment}
\concept{attachment function}
