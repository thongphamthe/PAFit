\name{plot.Full_PAFit_result}
\alias{plot.Full_PAFit_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting the estimated attachment function and node fitness
}
\description{
This function plots the estimated attachment function \eqn{A_k} and node fitness \eqn{eta_i}, together with additional information such as their confidence intervals or the estimated attachment exponent (\eqn{\alpha} when assuming \eqn{A_k = k^\alpha}). 
}
\usage{
\method{plot}{Full_PAFit_result}(x,
     net_stat                 ,
     true_f         = NULL    , plot             = "A"              , plot_bin   = TRUE ,
     line           = FALSE   , confidence       = TRUE             , high_deg_A = 1    ,
     high_deg_f     = 5       ,
     shade_point    = 0.5     , col_point        = "grey25"         , pch        = 16   ,
     shade_interval = 0.5     , col_interval     = "lightsteelblue" , label_x    = NULL , 
     label_y        = NULL    ,
     max_A          = NULL    , min_A            = NULL             , f_min      = NULL , 
     f_max          = NULL    , plot_true_degree = FALSE , 
     ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{Full_PAFit_result}, containing the estimated results from \code{\link{only_A_estimate}}, \code{\link{only_F_estimate}} or \code{\link{joint_estimate}}.
}
  \item{net_stat}{
An object of class \code{PAFit_data}, containing the summerized statistics.
}
  \item{true_f}{
Vector. Optional parameter for the true value of node fitnesses (only available in simulated datasets). If this parameter is specified and \code{plot == "true_f"}, a plot of estimated \eqn{\eta} versus true \eqn{\eta} is produced (after a suitable rescaling of the estimated \eqn{f}).  
}
  \item{plot}{
String. Indicates which plot is produced. 
\itemize{
\item If \code{"A"} then PA function is plotted. 
\item If \code{"f"} then the histogram of estimated fitness is plotted. 
\item If \code{"true_f"} then estimated fitness and true fitness are plotted together (require supplement of true fitness). 
}
Default value is \code{"A"}.
}
\item{plot_bin}{Logical. If \code{TRUE} then only the center of each bin is plotted. Default is \code{TRUE}.}
\item{line}{
Logical. Indicates whether to plot the line fitted from the log-linear model or not. Default value is \eqn{TRUE}.
}
\item{confidence}{
Logical. Indicates whether to plot the confidence intervals of \eqn{A_k} and \eqn{eta_i} or not. If confidence == TRUE, a 2-sigma confidence interval will be plotted at each \eqn{A_k} and \eqn{eta_i}.
}
\item{high_deg_A}{Integer. The estimated PA function is plotted starting from \code{high_deg_A}. Default value is \code{1}.}

\item{high_deg_f}{Integer. If \code{plot == "true_f"}, only nodes whose number of edges acquired is not less than \code{high_deg_f} are plotted. Default value is \code{5}.}

\item{col_point}{String. The name of the color of the points. Default value is \code{"black"}.}

\item{shade_point}{
Numeric. Value between 0 and 1. This is the transparency level of the points. Default value is \code{0.5}.
}
\item{pch}{Numeric. The plot symbol. Default value is \code{16}.}
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
#plot A
plot(result , net_stats , plot = "A")
true_A     <- c(1,result$estimate_result$center_k[-1])
lines(result$estimate_result$center_k + 1 , true_A , col = "red") # true line
legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
#plot true_f
plot(result, net_stats , net$fitness, plot = "true_f")
}
}
