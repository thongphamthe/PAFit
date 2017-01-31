\name{CreateDataCV}
\alias{CreateDataCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating cross validation data
}
\description{
A function to create cross-validation data.
}
\usage{
CreateDataCV(net                   , p          = 0.75 , G           = 50  , 
             net_type = "directed" , deg_thresh = 0    , exclude_end = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
A three-column matrix whose each row contains information of one edge in the form \code{(from_node id, to_node id, time_stamp)}. \code{from_node id} is the id of the source node. \code{to_node id} is the id of the destination node. \code{time_stamp} is the arrival time of the edge. \code{from_node id} and \code{to_node id} are assumed to be integers starting from \eqn{0}. \code{time_stamp} can be either numeric or string. The value of a time-stamp can be arbitrary, but we assume that a smaller time_stamp (regarded so by the \code{sort} function in \code{R}) represents an earlier arrival time.}
\item{p}{
Numeric between \code{0} and \code{1}. Indicates the ratio of number of new edges in the learning data to that of the full data. Default is \code{p = 0.75}.
}
\item{G}{
Integer. Number of bins. Default value is \\code{50}.
}
 \item{net_type}{
String. The type of the network: \code{"directed"} or \code{"undirected"}. Default is \code{"directed"}. 
}
\item{deg_thresh}{
Integer. We only consider nodes with number of acquired new edges at least this threshold. Default value is \code{0},i.e. all the nodes.
}
\item{exclude_end}{
Logical. If \code{TRUE}, then for the testing data, at each time-step we only consider the new edges that connect to nodes with the current degrees less than \eqn{deg\_max}, which is the maximum degree in the learning data. The motivation for this option is that in the learning phase, we can only learn the PA function up to \eqn{deg_max}, so it makes sense to limit the degree in the testing phase to \eqn{deg\_max}. From our experiences, this option does not matter. Default value is \code{FALSE}}.
}

\value{An object of class \code{"CV_Data"} containing the data needed for cross validation.
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
net      <- GenerateNet(N = 100 , m = 1 , mode = 1 , alpha = 1 , shape = 5 , rate = 5)
data_cv  <- CreateDataCV(net$graph)
summary(data_cv)
}

\keyword{fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{Bianconi-Barabasi model}
\keyword{Barabasi-Albert model}% __ONLY ONE__ keyword per line
