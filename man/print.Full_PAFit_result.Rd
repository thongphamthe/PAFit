\name{print.Full_PAFit_result}
\alias{print.Full_PAFit_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  printing information on the estimation result
}
\description{
  This function outputs simple information of the estimation result.
}
\usage{
  \method{print}{Full_PAFit_result}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class \code{Full_PAFit_result}, containing the estimated results from \code{\link{only_A_estimate}}, \code{\link{only_F_estimate}} or \code{\link{joint_estimate}}.
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
    print(result)
  }
}
