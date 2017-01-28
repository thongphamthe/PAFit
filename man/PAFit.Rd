\name{PAFit}
\alias{PAFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to estimate the attachment function and node fitness in temporal complex networks
}
\description{
From a \code{PAFit_data} object, which contains summary statistics summerized from the dataset, PAFit estimates the attachment function \eqn{A_k} and node fitness \eqn{\eta_i} by penalized log-likelihood maximization.  It also calculates confidence intervals for \eqn{A_k} and \eqn{\eta_i}.  Estimation of the attachment function in isolation (while fixing \eqn{\eta_i = 1}) can be specified by setting only_PA = TRUE.  Estimation of node fitness in isolation (fixing either \eqn{A_k = k} or \eqn{A_k = 1}) can be specified by setting \code{only_f = TRUE}.
}
\usage{
PAFit (net_stat, 
       only_PA        = FALSE       , only_f      = FALSE, 
       mode_f         = "Linear_PA" ,
       true_A         = NULL        , true_f      = NULL , s           = 10   , 
       auto_lambda    = TRUE        , r           = 0.01 , lambda      = 1    , 
       weight_PA_mode = 1           , auto_stop   = TRUE , stop_cond   = 10^-6, 
       iteration      = 200         , max_iter    = 2e+05, debug       = FALSE, 
       step_size      = 0.5         , q           = 1    , alpha_start = 1    , 
       normalized_f   = FALSE       , interpolate = TRUE , mode_reg_A  = 0    , 
       start_mode_A   = "Log_linear", 
       start_mode_f   = "Constant"  ,
       ...)
}
\arguments{
  \item{net_stat}{
An object of class "PAFit_data" containing all the summary statistics summerized from the data by the function GetStatistics
}
\item{only_PA}{
Logical. TRUE means that the attachment function \eqn{A_k} is estimated in isolation(fixing \eqn{\eta_i = 1}). Default is FALSE.
}
\item{only_f}{
Logical. TRUE means that the fitness function is estimated in isolation. Default is FALSE.
}
\item{mode_f}{
String. Only effective when only_f == TRUE. If mode_f == "Linear_PA" then \eqn{A_k = k} for \eqn{k \ge 1} and \eqn{A_0 = C}. If mode_f == "Constant_PA" then \eqn{A_k = 1} for all \eqn{k}. If mode_f == "Log_linear" then \eqn{A_k = k^\alpha} for \eqn{k \ge 1} and \eqn{A_0 = C}. The values of \eqn{\alpha} and \eqn{C} are estimated by MLE. Default values is "Linear_PA".
}
\item{true_A}{
Numeric vector. If the true_A is supplemented, then only fitness is estimated.
}
\item{true_f}{
Numeric vector. If the true_f is supplemented, then only PA is estimated.
}
\item{s}{
Numeric. The regularization parameter \eqn{s} for node fitness. Default value is \eqn{1}.
}

\item{auto_lambda}{
Logical. If auto_lambda == TRUE, lambda will be determined automatically from the data. Default is TRUE.
}
\item{r}{
Numeric. The regularization parameter \eqn{r} for the PA function. Default value is \eqn{0.01}.
}

\item{lambda}{
Numeric. The strength of the regularization for PA function. Ignored when auto_lambda == TRUE. Default value is \eqn{1}. lambda == 0 means no regularization for \eqn{A}. 
}
\item{weight_PA_mode}{
Integer. Indicates how the regularization terms for \eqn{A_k} are weighted. If weight_PA_mode == 0, the regularization term for \eqn{A_k} is weighted by the total number of edges connected to degree \eqn{k} nodes. If weight_PA_mode == 1, the regularization terms have uniform weights. Default value is \eqn{0}.
}
\item{auto_stop}{Logical. Indicates whether the algorithm stop automatically or not. Default is TRUE}
\item{stop_cond}{Numeric. If \eqn{auto_stop = TRUE}, the iterative algorithm stops when \eqn{abs(h(ii) - h(ii + 1)) / (abs(h(ii)) + 1) < stop_cond} where \eqn{h(ii)} is the value of the objective function (posterior probability in log-scale) at iteration \eqn{ii}. We recommend to choose stop_cond at most equal to \eqn{10^(- number of digits of h - 2)}, in order to ensure that when the algorithm stops, the increase in posterior probability is less than 1\% of the current posterior probability. Default is \eqn{10^-9}.}
\item{iteration}{
Numeric. The number of iterations. Ignored if auto_stop == TRUE. Default value is \eqn{20}. 
}
\item{max_iter}{
Numeric. The maximum number of iterations. Regardless of other settings, the algorithm will stop once the number of  iterations reaches this threshold. Default value is \eqn{2e+05}.
}
\item{debug}{Logical. if debug == TRUE, the value of the objective function \eqn{h} is printed out at each step. Defaule is FALSE.}
\item{step_size}{Numeric. A number between (0,1] to indicate the step-size of the quasi-Newton speedup. Ignored (no quasi-Newton speedup) if \eqn{q <= 1}. Defaule is 1.}
\item{q}{Integer. Indicates numbers of previous steps using in the quasi-Newton speedup. Ignored if \eqn{q <= 1}. Defaule is 0.5.}
\item{alpha_start}{Numeric. The starting value for alpha when we use the model \eqn{k^\alpha}. Default value is \eqn{1}.}
\item{normalized_f}{Logical. Indicates whether we should normalize the estimated value of \eqn{f} after estimation. Default value is FALSE.}
\item{interpolate}{Logical. Indicates whether we should perform interpolation for the missing values of the estimated \eqn{A_k}. The interpolation, if performed, is a linear regression on log-scale. Default value is TRUE.}
\item{mode_reg_A}{Integer. Indicate which regularization term is used for the PA function. For the regularization function used in the PLOS ONE and SR paper, use 0. Default value is \eqn{0}.}
\item{start_mode_A}{String. Takes one of two values: "Log_linear" (the initial PA function set to k^alpha_start) or "Random" (the initial function is randomly sampled from a uniform distribution). Default value is "Log_linear".}
\item{start_mode_f}{String. Takes one of two values: "Constant" (the initial node fitnesses are all set to 1) or "Random" (the initial node fitnesses are randomly sampled from a gamma distribution). Default value is "Constant".}


\item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\value{
an object of class "PAFit_result", which is a list. Some important fields of this object:
\item{A }{The estimated attachment function}
\item{k}{The corresponding degree}
\item{var_A}{Variances of the estimated A}
\item{linear_fit}{Result of fitting the log-linear model \eqn{log A_k = \alpha log k} to the estimated \eqn{A_k}}
\item{alpha}{The estimated attachment exponent of the log-linear model \eqn{A_k =k^\alpha}}
\item{weight_of_A}{The number of A in each bin}
\item{var_logA}{Variances of log(A)}        
\item{upper_A}{The upper value of the 2-sigma confidence interval of A}
\item{lower_A}{The lower value of the 2-sigma confidence interval of A}
\item{center_k}{The logarithmic center of the bins}
\item{theta}{Attachment value of the bins (before converting back to \eqn{A_k})}
\item{upper_bin}{The upper value of the 2-sigma confidence interval of theta}
\item{lower_bin}{The lower value of the 2-sigma confidence interval of theta}
\item{f}{The estimated node fitnesses \eqn{\eta}}
\item{var_f}{Variances of the estimated node fitnesses}
\item{upper_f}{The upper value of the 2-sigma confidence interval of \eqn{\eta}}
\item{lower_f}{The lower value of the 2-sigma confidence interval of \eqn{\eta}}
\item{objective_value}{Values of the objective function \eqn{h} (posterior probability in log-scale) recorded at each iteration}
}

\author{
Thong Pham
\email{thongpham@thongpham.net}
}

\references{
1. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Nonparametric Estimation of the Preferential Attachment Function in Complex Networks: Evidence of Deviations from Log Linearity, Proceedings of ECCS 2014, 141-153 (Springer International Publishing) (\url{http://dx.doi.org/10.1007/978-3-319-29228-1_13}).

2. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).

3. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{www.nature.com/articles/srep32558}).
}


\examples{
library("PAFit")
net        <- GenerateNet(N = 50,m = 1, mode = 1, alpha = 1, shape = 10, rate = 10)
net_stats  <- GetStatistics(net$graph)
result     <- PAFit(net_stats)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package }
\keyword{fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{Barabasi-Albert model}
\keyword{MM algorithims}
\keyword{Bianconi-Barabasi model}
\keyword{fitness}
\keyword{maximum likelihood function}% __ONLY ONE__ keyword per line
