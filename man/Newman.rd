\name{Newman}
\alias{Newman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Corrected Newman's method for estimating the preferential attachment function
}
\description{
This function implements a correction proposed in [1] of the original Newman's method in [2] to estimate the preferential attachment function. 
}
\usage{
  Newman(net_object                              , 
         net_stat    = get_statistics(net_object), 
         start       = 1                         , 
         interpolate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{net_object}{
    an object of class \code{PAFit_net} that contains the network.
}
  \item{net_stat}{
    An object of class \code{PAFit_data} which contains summerized statistics needed in estimation. This object is created by the function \code{\link{get_statistics}}. Default value is \code{ get_statistics(net_object)}.
  }
  \item{start}{Positive integer. The starting time from which the method is applied. Default value is \eqn{1}.}
  \item{interpolate}{
    Logical. If \code{TRUE} then all the gaps in the estimated PA function are interpolated by linear interpolating in logarithm scale. Default value is \code{FALSE}.
  }
}
\value{
  Outputs an \code{PA_result} object which contains the estimated attachment function. In particular, it contains the following field:
   \itemize{
    \item \code{k} and \code{A}: a degree vector and the estimated PA function.
    
    \item \code{center_k} and \code{theta}: when we perform binning, these are the centers of the bins and the estimated PA values for those bins. 
    \item \code{g}: the number of bins used.
    \item \code{alpha} and \code{ci}: \code{alpha} is the estimated attachment exponenet \eqn{\alpha} (when assume \eqn{A_k = k^\alpha}), while \code{ci} is the mean plus/minus two-standard-deviation interval.
    \item \code{loglinear_fit}: this is the fitting result when we estimate \eqn{\alpha}. 
}
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).
  
  2. Newman, M.. Clustering and preferential attachment in growing networks. Physical Review E. 2001;64(2):025102 (\url{https://journals.aps.org/pre/abstract/10.1103/PhysRevE.64.025102}).
}
\seealso{

See \code{\link{get_statistics}} for how to create summerized statistics needed in this function.
 
See \code{\link{Jeong}}, \code{\link{only_A_estimate}} for other methods to estimate the attachment function in isolation.
}
\examples{
  library("PAFit")
  net        <- generate_net(N = 1000 , m = 1 , mode = 1 , alpha = 1 , s = 0)
  net_stats  <- get_statistics(net)
  result     <- Newman(net, net_stats)
  summary(result)
  # true function
  true_A     <- result$center_k
  #plot the estimated attachment function
  plot(result , net_stats)
  lines(result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
}

\concept{preferential attachment}
\concept{attachment function}
