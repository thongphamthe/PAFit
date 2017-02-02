\name{PAFit}
\alias{PAFit}
\alias{PAFit_result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint inference of preferential attachment and node fitness by Minorize-Maximization algorithms
}
\description{
From a \code{PAFit_data} object, which contains summary statistics of the dataset, PAFit estimates the attachment function \eqn{A_k} and node fitness \eqn{\eta_i} by penalized log-likelihood maximization. It also infers the remaining uncertainties in the estimated results by approximating the confidence intervals of \eqn{A_k} and \eqn{\eta_i}.  

Estimation of either the attachment function or node fitness in isolation are also supported. Estimation of the PA function with \eqn{\eta_i = 1} can be specified by setting \code{only_PA = TRUE}.  Estimation of node fitness with either \eqn{A_k = k} or \eqn{A_k = 1} can be specified by setting \code{only_f = TRUE}.
}
\usage{
PAFit (net_stat, 
       only_PA        = FALSE       , only_f         = FALSE       , 
       mode_f         = "Linear_PA" ,
       true_A         = NULL        , true_f         = NULL        , 
       
       mode_reg_A     = 0           , weight_PA_mode = 1           ,
       s              = 10          , lambda         = 1           , 
       auto_lambda    = TRUE        ,  r             = 0.01        , 
       
       alpha_start    = 1           , start_mode_A   = "Log_linear", 
       start_mode_f   = "Constant"  ,
       
       auto_stop      = TRUE        , stop_cond      = 10^-7       , 
       iteration      = 200         , max_iter       = 2e+05       , 
       debug          = FALSE       , q              = 1           ,   
       step_size      = 0.5         ,
      
       normalized_f   = FALSE       , interpolate    = FALSE)
}
\arguments{
The parameters can be divided into five groups based on what they specify.

First group specifies basic instructions for the algorithm: 

  \item{net_stat}{
An object of class \code{PAFit_data} containing all the summary statistics summerized from the data by the function \code{\link{GetStatistics}}.
}
\item{only_PA}{
Logical. \code{TRUE} means that the attachment function \eqn{A_k} is estimated in isolation(fixing \eqn{\eta_i = 1}). Default is \code{FALSE}.
}
\item{only_f}{
Logical. \code{TRUE} means that the fitness function is estimated in isolation. Default is \code{FALSE}.
}
\item{mode_f}{
String. Possible values: \code{"Linear_PA"}, \code{"Constant_PA"} or \code{"Log_linear"}. In the first two cases, the PA function is fixed. If \code{mode_f == "Linear_PA"} then \eqn{A_k = k} for \eqn{k \ge 1} and \eqn{A_0 = 1}. If \code{mode_f == "Constant_PA"} then \eqn{A_k = 1} for all \eqn{k}. In the final case of \code{mode_f == "Log_linear"}, we set \eqn{A_k = k^\alpha} for \eqn{k \ge 1} and \eqn{A_0 = 1}. The value of \eqn{\alpha} is also estimated. Default values is \code{"Linear_PA"}.
}
\item{true_A}{
Numeric vector. User-supplemented value of the PA function. If \code{true_A} is supplemented, then only node fitnesses are estimated.
}
\item{true_f}{
Numeric vector. User-supplemented value of node fitnesses. If \code{true_f} is supplemented, then only the PA function is estimated.
}

Second group specifies the objective function, e.g. the weighting of PA, the regularization terms for PA and fitness: 

\item{mode_reg_A}{Integer. Possible values: \code{0}, \code{1} or \code{2}. Indicates which regularization term is used for the PA function. For the regularization function used in the PLOS ONE and SR paper, use \code{0}. Default value is \code{0}.}

\item{weight_PA_mode}{
Binary. Indicates how the regularization terms for \eqn{A_k} are weighted. If \code{weight_PA_mode == 0}, the regularization term for \eqn{A_k} is weighted by the total number of edges connected to degree \eqn{k} nodes. If \code{weight_PA_mode == 1}, the regularization terms have uniform weights. Default value is \eqn{0}.
}
\item{s}{
Positive numeric. The regularization parameter \code{s} for node fitness. Default value is \code{10}. In real-world situations, since we do not know beforehand the suitable value for \code{s}, we should choose its value by cross-validation. This is done in the \code{\link{performCV}} function.
}
\item{lambda}{
Non-negative numeric. The absolute strength of the regularization for PA function. Ignored when  \code{auto_lambda == TRUE}. Default value is  \code{1}.  \code{lambda == 0} means no regularization for PA. 
}
\item{auto_lambda}{
Logical. If \code{auto_lambda == TRUE}, \code{lambda} will be determined automatically from the data by \code{r}. Default is \code{TRUE}.
}
\item{r}{
Non-negative numeric. The regularization parameter \code{r} for the PA function indicates the relative strength of the regularization term. From \code{r}, the value of \code{lambda} is automatically determined if \code{auto_lambda == TRUE}. Default value is \code{0.01}. In real-world situations, since we do not know beforehand the suitable value for \code{r}, we should choose its value by cross-validation. This is done in the \code{\link{performCV}} function.
}



Third group specifies the initial value of PA and fitness: 

\item{alpha_start}{Non-negative numeric. The starting value for \eqn{\alpha} when we use the model \eqn{k^\alpha}. Default value is \code{1}.}
\item{start_mode_A}{String. Takes one of two values: \code{"Log_linear"} (the initial PA function set to k^alpha_start) or \code{"Random"} (the initial function is randomly sampled from a uniform distribution). Default value is \code{"Log_linear"}.}
\item{start_mode_f}{String. Takes one of two values: \code{"Constant"} (the initial node fitnesses are all set to \code{1}) or \code{"Random"} (the initial node fitnesses are randomly sampled from a gamma distribution). Default value is \code{"Constant"}.}

Fourth group concerns the iterative process:

\item{auto_stop}{Logical. Indicates whether the algorithm stop automatically or not. Default is \code{TRUE}}
\item{stop_cond}{Numeric. If \code{auto_stop = TRUE}, the iterative algorithm stops when \eqn{abs(h(ii) - h(ii + 1)) / (abs(h(ii)) + 1) < stop_cond} where \eqn{h(ii)} is the value of the objective function at iteration \eqn{ii}. We recommend to choose stop_cond at most equal to \eqn{10^(- number of digits of h - 2)}, in order to ensure that when the algorithm stops, the increase in posterior probability is less than 1\% of the current posterior probability. Default is \code{10^-7}.}
\item{iteration}{
Integer. The number of iterations. Ignored if \code{auto_stop == TRUE}. Default value is \code{20}. 
}
\item{max_iter}{
Integer. The maximum number of iterations. Regardless of other settings, the algorithm will stop once the number of  iterations reaches this threshold. Default value is \code{2e+05}.
}
\item{debug}{Logical. if \code{debug == TRUE}, the value of the objective function \eqn{h} is printed out at each step. Defaule is \code{FALSE}.}
\item{q}{Integer. Indicates numbers of previous steps using in the quasi-Newton speedup. Ignored if \eqn{q <= 1}. Defaule is \code{1}.}
\item{step_size}{Numeric. A number between \eqn{(0,1]} to indicate the step-size of the quasi-Newton speedup. Ignored (no quasi-Newton speedup) if \eqn{q <= 1}. Defaule is \code{0.5}.}

Final group gives some additional instructions: 

\item{normalized_f}{Logical. Indicates whether we should normalize the estimated value of node fitness after estimation. Default value is \code{FALSE}.}
\item{interpolate}{Logical. Indicates whether we should perform interpolation for the missing gaps in the estimated \eqn{A_k}. The interpolation, if performed, is a linear regression on log-scale. Default value is \code{FALSE}.}

}
\value{
an object of class \code{PAFit_result}, which is a list. Some important fields can be divided into five groups. 


The first group gives the estimated preferential attachment function:


\item{k}{The observed degree vector}
\item{A }{The estimated attachment function corresponding to \code{k}}
\item{center_k}{The logarithmic center of the bins}
\item{theta}{Preferential attachment value corresponding to \code{center_k} (before mapping back to \eqn{A_k})}
\item{weight_of_A}{The number of \code{A} in each bin}
\item{loglinear_fit}{Result of fitting the log-linear model \eqn{log A_k = \alpha log k + C} to the estimated \eqn{A_k}}
\item{alpha}{The estimated attachment exponent of the log-linear model \eqn{A_k =k^\alpha}}
\item{ci}{The confidence interval of the attachment exponent. It is two-sigma. When \code{mode_f != "Log_linear"}, this confidence interval is estimated from the log_linear fit (fitting \eqn{log k} to \eqn{log A_k}) using \code{confint} function, so it has a popular meaning as a \eqn{95}-percentage confidence interval.}
\item{alpha_series}{The series of \eqn{\alpha} over iterations if \code{mode_f == "Log_linear"}}     

The second group gives the confidence intervals of the estimated PA function:

\item{var_A}{Variances of the estimated A}
\item{var_logA}{Variances of \eqn{log(A)}}   
\item{upper_A}{The upper value of the two-sigma confidence interval of A}
\item{lower_A}{The lower value of the two-sigma confidence interval of A}
\item{upper_bin}{The upper value of the two-sigma confidence interval of theta}
\item{lower_bin}{The lower value of the two-sigma confidence interval of theta}


The third group gives the estimated node fitnesses: 

\item{f}{The estimated node fitnesses \eqn{\eta}}


The fourth group gives the confidence intervals of the estimated node fitnesses:

\item{var_f}{Variances of the estimated node fitnesses}
\item{upper_f}{The upper value of the two-sigma confidence interval of node fitness \eqn{\eta}}
\item{lower_f}{The lower value of the two-sigma confidence interval of node fitness \eqn{\eta}}



The final group gives additional information on the iterative process:

\item{objective_value}{Values of the objective function \eqn{h} recorded at each iteration}


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

\seealso{\code{\link{Newman_corrected}}, \code{\link{Jeong}}, \code{\link{CreateDataCV}}, \code{\link{performCV}}}

\examples{
library("PAFit")
net        <- GenerateNet(N = 50,m = 10, mode = 1, alpha = 0.5, shape = 100, rate = 100)
net_stats  <- GetStatistics(net$graph)
result     <- PAFit(net_stats, r = 0.01, s = 100)
summary(result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package }
\keyword{fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{\enc{Barabási}{Barabasi}-Albert model}
\keyword{Bianconi-\enc{Barabási}{Barabasi} model}
\keyword{fitness}
\keyword{maximum likelihood function}% __ONLY ONE__ keyword per line
