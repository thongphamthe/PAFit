\name{plot.PA_result}
\alias{plot.PA_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Plotting the preferential attachment function
}
\description{
    This function plots the estimated attachment function from the corrected Newman's method or the Jeong's method. Its also plots additional information such as the estimated attachment exponenent (\eqn{\alpha} when assuming \eqn{A_k = k^\alpha}). 
}
\usage{
\method{plot}{PA_result}(x, 
     net_stat,
     plot_bin    = TRUE   ,
     high_deg    = NULL   ,  
     line        = FALSE  , 
     col_point   = "black",
     shade_point = 0.5    , 
     max_A       = NULL   , 
     min_A       = NULL   , 
     label_x     = NULL   , 
     label_y     = NULL   ,
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class \code{PA_result}, containing the estimated PA function and the estimated attachment exponenet from either \code{\link{Newman}} or \code{\link{Jeong}} functions. 
  }
  \item{net_stat}{
    An object of class \code{PA_data}, containing the summerized statistics. This object is created from the function \code{\link{GetStatistics}}.
  }
  \item{plot_bin}{Logical. If \code{TRUE} then only the center of each bin is plotted. Default is \code{TRUE}.}
  \item{high_deg}{Integer. Specifies the starting degree from which \eqn{A_k} is plotted. If this parameter is specified, the estimated PA function is plotted from \code{k = high_deg}}
  \item{line}{
    Logical. Indicates whether to plot the line fitted from the log-linear model or not. Default value is \code{FALSE}.
  }
  \item{col_point}{String. The name of the color of the points. Default value is \eqn{"black"}.}
  \item{shade_point}{
    Numeric. Value between \code{0} and \code{1}. This is the transparency level of the points. Default value is \code{0.5}.
  }
  \item{max_A}{Numeric. Specify the maximum of the axis of PA.}
  \item{min_A}{Numeric. Specify the minimum of the axis of PA.}
  \item{label_x}{String. The label of x-axis.}
  \item{label_y}{String. The label of y-axis.}
  \item{\dots}{
   Additional parameters to pass onto the \code{plot} function.
  }
}
\value{
  Outputs the desired plot.
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).
}

\examples{
  library("PAFit")
  net        <- GenerateNet(N = 1000 , m = 1 , mode = 1 , alpha = 1 , shape = 0)
  net_stats  <- GetStatistics(net$graph)
  result     <- Newman(net_stats)
  # true function
  true_A     <- result$center_k
  #plot the estimated PA function
  plot(result , net_stats)
  lines(result$center_k + 1, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
}
