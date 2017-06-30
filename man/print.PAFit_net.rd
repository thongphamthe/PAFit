\name{print.PAFit_net}
\alias{print.PAFit_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Printing simple information of a \code{PAFit_net} object
}
\description{
  This function outputs simple information of a \code{PAFit_net} object.
}
\usage{
  \method{print}{PAFit_net}(x,
                            ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class \code{PAFit_net}.
  }
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}


\value{
  Outputs simple information of the network. 
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}


\examples{
  library("PAFit")
  # a network from Bianconi-Barabasi model
  net        <- generate_BB(N = 50 , m = 10 , s = 10)
  print(net)
}
