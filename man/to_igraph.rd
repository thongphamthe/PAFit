\name{to_igraph}
\alias{to_igraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Convert a PAFit_net object to an igraph object
}
\description{
  This function converts a \code{PAFit_net} object to an \code{igraph} object (of package \pkg{igraph}).
}
\usage{
to_igraph(net_object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{net_object}{
An object of class \code{PAFit_net}.
}
}

\value{
The function returns an \code{igraph} object.
}

\author{
Thong Pham \email{thongpham@thongpham.net}
}


\examples{
library("PAFit")
# a network from Bianconi-Barabasi model
net          <- generate_BB(N = 50 , m = 10 , s = 10)
igraph_graph <- to_igraph(net)
}
