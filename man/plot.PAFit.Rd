\name{plot.PAFit}
\alias{plot.PAFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot the estimated attachment function and node fitness
}
\description{
This function plots the estimated attachment function \eqn{A_k} and node fitness \eqn{eta_i}, together with additional information such as their confidence intervals. 
}
\usage{
\method{plot}{PAFit}(x,data,true_f = NULL, plot = c("A","f","true_f"), plot_bin = TRUE,
         line = FALSE, confidence = TRUE, high_deg = NULL, 
         shade_point = 0.5, shade_interval = 0.5, 
         max_A = NULL, 
         min_A = NULL, f_min = NULL, 
         f_max = NULL, plot_true_degree = FALSE, 
         label_x = NULL, label_y = NULL,
         col_interval = "lightsteelblue",
         col_point = "black",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "PAFit", containing the result
}
  \item{data}{
An object of class "PAFitData", containing the summerized statistics.
}
  \item{true_f}{
Vector. Optional parameter for the true value of node fitnesses (only available in simulated datasets). If this parameter is specified and plot == "true_f", a plot of estimated \eqn{\eta} versus true \eqn{\eta} is produced (after a suitable rescaling of the estimated \eqn{f}).  
}
  \item{plot}{
String. Indicates which plot is produced.if \eqn{"A"} then PA function is plotted. If \eqn{"f"} then estimated fitness is plotted. If \eqn{"true_f"} then estimated fitness and true fitness are plotted together (require supplement of true fitness). Default value is \eqn{"A"}.
}
\item{plot_bin}{Logical. If TRUE then only the center of each bin is plotted. Default is \eqn{TRUE}.}
\item{line}{
Logical. Indicates whether to plot the line fitted from the log-linear model or not. Default value is \eqn{TRUE}.
}
\item{confidence}{
Logical. Indicates whether to plot the confidence intervals of \eqn{A_k} and \eqn{eta_i} or not. If confidence == TRUE, a 2-sigma confidence interval will be plotted at each \eqn{A_k} and \eqn{eta_i}.
}
\item{high_deg}{Integer. If this parameter is specified, only nodes whose number of edges acquired is not less than \eqn{high_deg} is plotted.}
\item{shade_point}{
Numeric. Value between 0 and 1. This is the transparency level of the points. Default value is \eqn{0.5}.
}
\item{shade_interval}{
Numeric. Value between 0 and 1. This is the transparency level of the confidence intervals. Default value is \eqn{0.5}.
}
\item{max_A}{Numeric. Specify the maximum of the axis of PA.}
\item{min_A}{Numeric. Specify the minimum of the axis of PA.}
\item{f_min}{Numeric. Specify the minimum of the axis of fitness.}
\item{f_max}{Numeric. Specify the maximum of the axis of fitness.}
\item{plot_true_degree}{Logical. The degree of each node is plotted or not.}
\item{label_x}{String. The label of x-axis.}
\item{label_y}{String. The label of y-axis.}
\item{col_interval}{String. The name of the color of the confidence intervals. Default value is \eqn{"lightsteelblue"}.}
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
data   <- GenerateNet(N = 100,m = 1,mode = 1, alpha = 1, shape = 5, rate = 5)
stats  <- GetStatistics(data$graph,deg_thresh = 1, Binning = TRUE, G = 50)
result <- PAFit(stats,stop_cond = 10^-3)
#plot A
plot(result,stats,plot = "A")
#plot f
plot(result,stats,plot = "f")
#plot true_f
plot(result,stats,data$fitness, plot = "true_f")
}
