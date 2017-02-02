\name{summary.PAFit_result}
\alias{summary.PAFit_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing a \code{PAFit_result} object
}
\description{
This function summarizes information of a \code{PAFit_result} object. 
}
\usage{
\method{summary}{PAFit_result}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{PAFit_result}, containing the estimation result of \code{PAFit}.
}
\item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\value{
Outputs summary information.
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
net        <- GenerateNet(N = 50,m = 10, mode = 1, alpha = 0.5, shape = 100, rate = 100)
net_stats  <- GetStatistics(net$graph)
result     <- PAFit(net_stats, r = 0.01, s = 100)
summary(result)
}
