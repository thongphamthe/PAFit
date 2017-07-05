\name{summary.PAFit_result}
\alias{summary.PAFit_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Output summary information on the estimation result stored in a \code{PAFit_result} object
}
\description{
  This function outputs summary information of the estimation result stored in a \code{PAFit_result} object. This object is stored in the field \code{$estimate_result} of a \code{Full_PAFit_result} object, which in turn is the returning value of \code{\link{only_A_estimate}}, \code{\link{only_F_estimate}} or \code{\link{joint_estimate}}. 
}
\usage{
  \method{summary}{PAFit_result}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An object of class \code{PAFit_result}.
  }
  
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}


\value{
  Outputs summary information on the estimation result.
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
    result     <- joint_estimate(net, net_stats)
    summary(result$estimate_result)
  }
}
