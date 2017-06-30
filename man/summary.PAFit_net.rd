\name{summary.PAFit_net}
\alias{summary.PAFit_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summary information of a \code{PAFit_net} object
}
\description{
  This function outputs summary information of a \code{PAFit_net} object.
}
\usage{
  \method{summary}{PAFit_net}(object,
                           ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An object of class \code{PAFit_net}.
  }
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}


\value{
  Outputs summary information of the network. 
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}


\examples{
  library("PAFit")
  # a network from Bianconi-Barabasi model
  net        <- generate_BB(N = 50 , m = 10 , s = 10)
  summary(net)
}
