\name{plot.PA_result}
\alias{plot.PA_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A function to plot the estimated attachment function from Newman's or Jeong's method. 
}
\description{
  This function plots the estimated attachment function \eqn{A_k} together with additional information such as the estimated attachment exponenent (estimated \eqn{\alpha} when assuming \eqn{A_k = k^\alpha}). 
}
\usage{
  \method{plot}{PA_result}(x, net_stat,
                           plot_bin = TRUE,
                           high_deg = NULL,  
                           line = FALSE, 
                           shade_point = 0.5, 
                           max_A = NULL, 
                           min_A = NULL, 
                           label_x = NULL, 
                           label_y = NULL,
                           col_point = "black",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class "PA_result", containing the result
  }
  \item{net_stat}{
    An object of class "PA_data", containing the summerized statistics.
  }
  \item{plot_bin}{Logical. If TRUE then only the center of each bin is plotted. Default is \eqn{TRUE}.}
  \item{high_deg}{Integer. If this parameter is specified, the estimated PA function is plotted from \eqn{k = high_deg}}
  \item{line}{
    Logical. Indicates whether to plot the line fitted from the log-linear model or not. Default value is \eqn{TRUE}.
  }
  \item{shade_point}{
    Numeric. Value between 0 and 1. This is the transparency level of the points. Default value is \eqn{0.5}.
  }
  \item{max_A}{Numeric. Specify the maximum of the axis of PA.}
  \item{min_A}{Numeric. Specify the minimum of the axis of PA.}
  \item{label_x}{String. The label of x-axis.}
  \item{label_y}{String. The label of y-axis.}
  \item{col_point}{String. The name of the color of the points. Default value is \eqn{"black"}.}
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
\references{
  1. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Nonparametric Estimation of the Preferential Attachment Function in Complex Networks: Evidence of Deviations from Log Linearity, Proceedings of ECCS 2014, 141-153 (Springer International Publishing) (\url{http://dx.doi.org/10.1007/978-3-319-29228-1_13}).
  
  2. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).
  
  3. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{www.nature.com/articles/srep32558}).
}

\examples{
  library("PAFit")
  net        <- GenerateNet(N = 1000 , m = 1 , mode = 1 , alpha = 1 , shape = 0)
  net_stats  <- GetStatistics(net$graph)
  result     <- Newman_corrected(net_stats)
  #plot A
  plot(result , net_stats)
}
