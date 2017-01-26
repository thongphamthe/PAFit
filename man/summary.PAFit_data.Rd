\name{summary.PAFit_data}
\alias{summary.PAFit_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to summarize a PAFit_data object 
}
\description{
This function summarizes information of a PAFit_data object. 
}
\usage{
\method{summary}{PAFit_data}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class "PAFit_data", containing the summary statistics needed in estimation.
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
data   <- GenerateNet(N = 100,m = 1,mode = 1, alpha = 1, shape = 5, rate = 5)
stats  <- GetStatistics(data$graph)
summary(stats)
}
