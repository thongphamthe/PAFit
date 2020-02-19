\name{plot.PA_result}
\alias{plot.PA_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Plotting the estimated attachment function
}
\description{
    This function plots the estimated attachment function from the corrected Newman's method or the Jeong's method. Its also plots additional information such as the estimated attachment exponenent (\eqn{\alpha} when assuming \eqn{A_k = k^\alpha}). 
}
\usage{
\method{plot}{PA_result}(x, 
     net_stat,
     plot_bin    = TRUE   ,
     high_deg    = 1      ,  
     line        = FALSE  , 
     col_point   = "black",
     shade_point = 0.5    , 
     pch         = 16     ,
     max_A       = NULL   , 
     min_A       = NULL   , 
     label_x     = NULL   , 
     label_y     = NULL   ,
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class \code{PA_result}, containing the estimated attachment function and the estimated attachment exponenet from either \code{\link{Newman}} or \code{\link{Jeong}} functions. 
  }
  \item{net_stat}{
    An object of class \code{PA_data}, containing the summerized statistics. This object is created from the function \code{\link{get_statistics}}.
  }
  \item{plot_bin}{Logical. If \code{TRUE} then only the center of each bin is plotted. Default is \code{TRUE}.}
  \item{high_deg}{Integer. Specifies the starting degree from which \eqn{A_k} is plotted. If this parameter is specified, the estimated attachment function is plotted from \code{k = high_deg}}
  \item{line}{
    Logical. Indicates whether to plot the line fitted from the log-linear model or not. Default value is \code{FALSE}.
  }
  \item{col_point}{String. The name of the color of the points. Default value is \eqn{"black"}.}
  \item{shade_point}{
    Numeric. Value between \code{0} and \code{1}. This is the transparency level of the points. Default value is \code{0.5}.
  }
\item{pch}{Numeric. The plot symbol. Default value is \code{16}.}  
  \item{max_A}{Numeric. Specify the maximum of the horizontal axis.}
  \item{min_A}{Numeric. Specify the minimum of the horizontal axis.}
  \item{label_x}{String. The label of x-axis. If \code{NULL}, then \code{"Degree k"} is used.}
  \item{label_y}{String. The label of y-axis. If \code{NULL}, then \code{"Attachment function"} is used.}
    \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}
\value{
  Outputs the desired plot.
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}

\examples{
  library("PAFit")
  net        <- generate_net(N = 1000 , m = 1 , mode = 1 , alpha = 1 , s = 0)
  net_stats  <- get_statistics(net)
  result     <- Newman(net, net_stats)
  # true function
  true_A     <- result$center_k
  # plot the estimated attachment function
  plot(result , net_stats)
  lines(result$center_k, true_A, col = "red") # true attachment function
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
}
