\name{plot.PAFit_result}
\alias{plot.PAFit_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting the estimated attachment function and node fitness
}
\description{
This function plots the estimated attachment function \eqn{A_k} and node fitness \eqn{eta_i}, together with additional information such as their confidence intervals or the estimated attachment exponent (\eqn{\alpha} when assuming \eqn{A_k = k^\alpha}). 
}
\usage{
\method{plot}{PAFit_result}(x,
     net_stat              ,
     true_f      = NULL    , plot             = "A"   , plot_bin     = TRUE             ,
     line        = FALSE   , confidence       = TRUE  , high_deg     = 1                , 
     shade_point = 0.5     , shade_interval   = 0.5   , col_interval = "lightsteelblue" ,
     col_point   = "black" , label_x          = NULL  , label_y      = NULL             ,
     max_A       = NULL    , min_A            = NULL  , f_min        = NULL             , 
     f_max       = NULL    , plot_true_degree = FALSE , 
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{PAFit_result}, containing the estimated results from \code{\link{PAFit}}.
}
  \item{net_stat}{
An object of class \code{PAFit_data}, containing the summerized statistics.
}
  \item{true_f}{
Vector. Optional parameter for the true value of node fitnesses (only available in simulated datasets). If this parameter is specified and \code{plot == "true_f"}, a plot of estimated \eqn{\eta} versus true \eqn{\eta} is produced (after a suitable rescaling of the estimated \eqn{f}).  
}
  \item{plot}{
String. Indicates which plot is produced. If \code{"A"} then PA function is plotted. If \code{"f"} then estimated fitness is plotted. If \code{"true_f"} then estimated fitness and true fitness are plotted together (require supplement of true fitness). Default value is \code{"A"}.
}
\item{plot_bin}{Logical. If TRUE then only the center of each bin is plotted. Default is \code{TRUE}.}
\item{line}{
Logical. Indicates whether to plot the line fitted from the log-linear model or not. Default value is \eqn{TRUE}.
}
\item{confidence}{
Logical. Indicates whether to plot the confidence intervals of \eqn{A_k} and \eqn{eta_i} or not. If confidence == TRUE, a 2-sigma confidence interval will be plotted at each \eqn{A_k} and \eqn{eta_i}.
}
\item{high_deg}{Integer. If \code{plot == "A"}, the estimated PA function is plotted starting from \code{high_deg}, and \eqn{A_{high_deg}} is normalized to \code{1}. If \code{plot == "f"} or \code{plot == "true_f"}, only nodes whose number of edges acquired is not less than \code{high_deg} are plotted. Default value is \code{1}.}

\item{shade_point}{
Numeric. Value between 0 and 1. This is the transparency level of the points. Default value is \code{0.5}.
}
\item{shade_interval}{
Numeric. Value between 0 and 1. This is the transparency level of the confidence intervals. Default value is \code{0.5}.
}
\item{max_A}{Numeric. Specify the maximum of the axis of PA.}
\item{min_A}{Numeric. Specify the minimum of the axis of PA.}
\item{f_min}{Numeric. Specify the minimum of the axis of fitness.}
\item{f_max}{Numeric. Specify the maximum of the axis of fitness.}
\item{plot_true_degree}{Logical. The degree of each node is plotted or not.}
\item{label_x}{String. The label of x-axis.}
\item{label_y}{String. The label of y-axis.}
\item{col_interval}{String. The name of the color of the confidence intervals. Default value is \code{"lightsteelblue"}.}
\item{col_point}{String. The name of the color of the points. Default value is \code{"black"}.}
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
net        <- GenerateNet(N = 50 , m = 1 , mode = 1 , alpha = 1 , shape = 10 , rate = 10)
net_stats  <- GetStatistics(net$graph)
result     <- PAFit(net_stats)
#plot A
plot(result , net_stats , plot = "A")
#plot f
plot(result , net_stats , plot = "f")
#plot true_f
plot(result , net_stats , net$fitness, plot = "true_f")
}
