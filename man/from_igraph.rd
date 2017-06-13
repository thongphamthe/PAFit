\name{from_igraph}
\alias{from_igraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert an igraph object to a PAFit_net object
}
\description{
  This function converts an \code{igraph} object (of package \pkg{igraph}) to a \code{PAFit_net} object.
}
\usage{
  from_igraph(net)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
    An object of class \code{igraph}.
  }
}

\value{
  The function returns a \code{PAFit_net} object.
}

\author{
  Thong Pham \email{thongpham@thongpham.net}
}


\examples{
  library("PAFit")
  # a network from Bianconi-Barabasi model
  net          <- generate_BB(N = 50 , m = 10 , s = 10)
  igraph_graph <- to_igraph(net)
  back         <- from_igraph(igraph_graph)
}
