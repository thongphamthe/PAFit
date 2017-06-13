\name{to_networkDynamic}
\alias{to_networkDynamic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert a PAFit_net object to a networkDynamic object
}
\description{
  This function converts a \code{PAFit_net} object to a \code{networkDynamic} object (of package \pkg{networkDynamic}).
}
\usage{
  to_networkDynamic(net_object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{net_object}{
    An object of class \code{PAFit_net}.
  }
}

\value{
  The function returns a \code{networkDynamic} object.
}

\author{
  Thong Pham \email{thongpham@thongpham.net}
}


\examples{
  library("PAFit")
  # a network from Bianconi-Barabasi model
  net          <- generate_BB(N = 50 , m = 10 , s = 10)
  nD_graph     <- to_networkDynamic(net)
}
