\name{print.PA_result}
\alias{print.PA_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Printing a \code{PA_result} object
}
\description{
  This function prints a \code{PA_result} object's information.
}
\usage{
\method{print}{PA_result}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{PA_result}, containing the estimated preferential attachment from either \code{\link{Newman_corrected}} or \code{\link{Jeong}} functions.
}
   \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}

\value{
Outputs some simple information. 
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
net        <- GenerateNet(N = 1000 , m = 30 , mode = 1 , alpha = 1 , shape = 0) # no fitness
net_stats  <- GetStatistics(net$graph)
result     <- Newman_corrected(net_stats)
print(result)
}
