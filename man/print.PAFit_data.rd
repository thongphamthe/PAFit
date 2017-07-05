\name{print.PAFit_data}
\alias{print.PAFit_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Printing simple information on the statistics of the network stored in a \code{PAFit_data} object
}
\description{
  This function prints simple information of the statistics stored in a \code{PAFit_data} object. This object is the returning value of \code{\link{get_statistics}}.
}
\usage{
  \method{print}{PAFit_data}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class \code{PAFit_data}.
  }
  
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}


\value{
  Prints simple information of the network statistics.
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
    print(net_stats)
  }
}
