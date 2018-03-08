\name{print.PA_result}
\alias{print.PA_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Printing information of the estimated attachment function
}
\description{
  This function outputs simple information of the estimated attachment function from the corrected Newman's method or the Jeong's method.
}
\usage{
  \method{print}{PA_result}(x, 
                              ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class \code{PA_result}, containing the estimated attachment function and the estimated attachment exponenet from either \code{\link{Newman}} or \code{\link{Jeong}} functions. 
  }
  \item{\dots}{
    Additional parameters to pass onto the \code{plot} function.
  }
}
\value{
  Simple information of the estimated attachment function.
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\examples{
  library("PAFit")
  net        <- generate_net(N = 1000 , m = 1 , mode = 1 , alpha = 1 , s = 0)
  net_stats  <- get_statistics(net)
  result     <- Newman(net, net_stats)
  print(result)
}
