\name{summary.PAFit_data}
\alias{summary.PAFit_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Output summary information on the statistics of the network stored in a \code{PAFit_data} object
}
\description{
  This function outputs summary information of the statistics stored in a \code{PAFit_data} object. This object is the returning value of \code{\link{get_statistics}}.
}
\usage{
  \method{summary}{PAFit_data}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An object of class \code{PAFit_data}.
  }
  
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}


\value{
  Outputs summary information of the network statistics.
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\examples{
  ## Since the runtime is long, we do not let this example run on CRAN
  \dontrun{
    library("PAFit")
    set.seed(1)
    # a network from Bianconi-Barabasi model
    net        <- generate_BB(N        = 1000 , m             = 50 , 
                              num_seed = 100  , multiple_node = 100,
                              s        = 10)
    net_stats  <- get_statistics(net)
    summary(net_stats)
  }
}
