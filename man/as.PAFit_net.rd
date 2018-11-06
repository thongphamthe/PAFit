\name{as.PAFit_net}
\alias{as.PAFit_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Converting an edgelist matrix to a PAFit_net object
}
\description{
  This function converts a graph stored in an edgelist matrix format to a \code{PAFit_net} object.
}
\usage{
as.PAFit_net(graph, type = "directed", PA = NULL, fitness = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{graph}{
An edgelist matrix. Each row is assumed to be of the form (\code{from_node_id} \code{to_node_id} \code{time_stamp}). For a directed network ,\code{from_node_id} is the id of the source node and \code{to_node_id} is the id of the destination node. For an undirected network, the order is ignored and \code{from_node_id} and \code{to_node_id} are the ids of two ends. \code{time_stamp} is the arrival time of the edge. \code{from_node_id} and \code{to_node_id} are assumed to be integers that are at least \eqn{0}. The whole ids need not to be contiguous.

To register a new node \eqn{i} at time \eqn{t} without any edge, add a row with format (\code{i -1 t}). This works for both undirected and directed networks.

\code{time_stamp} can be either numeric or string. The value of a time-stamp can be arbitrary, but we assume that a smaller time_stamp (regarded so by the \code{sort} function in \code{R}) represents an earlier arrival time. Examples of time-stamps that satisfy this assumption are the integer \code{0:T}, the string format `yyyy-mm-dd', and the POSIX time.

}
\item{type}{
String. Indicates whether the network is \code{"directed"} or \code{"undirected"}.
}
\item{PA}{
Numeric vector. Contains the PA function. Default value is \code{NULL}.
}
\item{fitness}{
Numeric vector. Contains node fitnesses. Default value is \code{NULL}.
}
}

\value{
An object of class \code{PAFit_net}
}

\author{
Thong Pham \email{thongpham@thongpham.net}
}


\examples{
library("PAFit")
# a network from Bianconi-Barabasi model
net        <- generate_BB(N = 50 , m = 10 , s = 10)
as.PAFit_net(net$graph)
}
