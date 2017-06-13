\name{from_networkDynamic}
\alias{from_networkDynamic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert a networkDynamic object to a PAFit_net object
}
\description{
  This function converts a \code{networkDynamic} object (of package \pkg{networkDynamic}) to a \code{PAFit_net} object.
}
\usage{
  from_networkDynamic(net)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net}{
    An object of class \code{networkDynamic}.
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
nD_graph     <- to_networkDynamic(net)
back         <- from_networkDynamic(nD_graph)
}
