\name{only_A_estimate}
\alias{only_A_estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimating the attachment function in isolation by PAFit method  
}
\description{
This function estimates the attachment function \eqn{A_k} by PAFit method. The method has a hyper-parameter \eqn{r}. It first performs a cross-validation step to select the optimal parameter \eqn{r} for the regularization of \eqn{A_k}, then uses that \eqn{r} to estimate the attachment function with the full data. 
}
\usage{
only_A_estimate(net_object                             , 
                net_stat   = get_statistics(net_object), 
                p          = 0.75                      ,
                stop_cond  = 10^-8                     , 
                mode_reg_A = 0                         ,
                MLE        = FALSE                     ,
               ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{net_object}{
    an object of class \code{PAFit_net} that contains the network.
  }
  \item{net_stat}{
    An object of class \code{PAFit_data} which contains summerized statistics needed in estimation. This object is created by the function \code{\link{get_statistics}}. The default value is \code{get_statistics(net_object)}.
  }
  \item{p}{Numeric. This is the ratio of the number of new edges in the learning data to that of the full data. The data is then divided into two parts: learning data and testing data based on \code{p}. The learning data is used to learn the node fitnesses and the testing data is then used in cross-validation. Default value is \code{0.75}.}
\item{stop_cond}{Numeric. The iterative algorithm stops when \eqn{abs(h(ii) - h(ii + 1)) / (abs(h(ii)) + 1) < stop.cond} where \eqn{h(ii)} is the value of the objective function at iteration \eqn{ii}. We recommend to choose \code{stop.cond} at most equal to \eqn{10^(- number of digits of h - 2)}, in order to ensure that when the algorithm stops, the increase in posterior probability is less than 1\% of the current posterior probability. Default is \code{10^-8}. This threshold is good enough for most applications.}

\item{mode_reg_A}{Binary. Indicates which regularization term is used for \eqn{A_k}:
\itemize{
\item \code{0}: This is the regularization term used in Ref. 1 and 2. Please refer to Eq. (4) in the tutorial for the definition of the term. It approximately enforces the power-law form \eqn{A_k = k^\alpha}. This is the default value. 
\item \code{1}: Unlike the default, this regularization term exactly enforces the functional form \eqn{A_k = k^\alpha}. Please refer to Eq. (6) in the tutorial for the definition of the term. Its main drawback is it is significantly slower to converge, while its gain over the default one is marginal in most cases.  
}
}
\item{MLE}{Logical. If \code{TRUE}, then not perform cross-validation and estimate the PA function with \code{r = 0}, i.e., maximum likelihood estimation. Default is \code{FALSE}. One might want to set this option to \code{TRUE} when one believes that there are sufficient data to get a reasonable MLE result, or when one wants to compare the default, regularized result with the MLE result.}
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}

\value{
   Outputs a \code{Full_PAFit_result} object, which is a list containing the following fields:
  \itemize{
  \item \code{cv_data}: a \code{CV_Data} object which contains the cross-validation data. This is the final  Normally the user does not need to pay attention to this data. \code{NULL} if \code{MLE = TRUE}.
  
 \item \code{cv_result}: a \code{CV_Result} object which contains the cross-validation result. Normally the user does not need to pay attention to this data. \code{NULL} if \code{MLE = TRUE}.
 
 \item \code{estimate_result}: this is a \code{PAFit_result} object which contains the estimated PA function and its confidence interval. It also includes the estimated attachment exponenent \eqn{\alpha} (assuming the model \eqn{A_k = k^\alpha}) in the field \code{alpha}, and the confidence interval of \eqn{\alpha} (in the field \code{ci}) when possible. In particular, the important fields are:
     \itemize{
     \item \code{ratio}: this is the selected value for the hyper-parameter \eqn{r}.
    \item \code{k} and \code{A}: a degree vector and the estimated PA function.
    \item \code{var_A}: the estimated variance of \eqn{A}.
    \item \code{var_logA}: the estimated variance of \eqn{log A}.
    \item \code{upper_A}: the upper value of the interval of two standard deviations around \eqn{A}.
    \item \code{lower_A}: the lower value of the interval of two standard deviations around \eqn{A}.
    
    \item \code{center_k} and \code{theta}: when we perform binning, these are the centers of the bins and the estimated PA values for those bins. \code{theta} is similar to \code{A} but with duplicated values removed.
     \item \code{var_bin}: the variance of \code{theta}. Same as \code{var_A} but with duplicated values removed.
    \item \code{upper_bin}: the upper value of the interval of two standard deviations around \code{theta}. Same as \code{upper_A} but with duplicated values removed.
    \item \code{lower_lower}: the lower value of the interval of two standard deviations around \code{theta}. Same as \code{lower_A} but with duplicated values removed.
    \item \code{g}: the number of bins used.
    \item \code{alpha} and \code{ci}: \code{alpha} is the estimated attachment exponenet \eqn{\alpha} (when assume \eqn{A_k = k^\alpha}), while \code{ci} is the confidence interval.
    \item \code{loglinear_fit}: this is the fitting result when we estimate \eqn{\alpha}. 
    \item \code{objective_value}: values of the objective function over iterations in the final run with the full data.
    \item \code{diverge_zero}: logical value indicates whether the algorithm diverged in the final run with the full data.
}
}
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).
  
  2. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{http://www.nature.com/articles/srep32558}).
}
\seealso{
  See \code{\link{get_statistics}} for how to create summerized statistics needed in this function.
  
  See \code{\link{Newman}} and \code{\link{Jeong}} for other methods to estimate the attachment function \eqn{A_k} in isolation.
  
}

\examples{
\dontrun{
  library("PAFit")
  set.seed(1)
  #### Example 1: Linear preferential attachment  #########
  # a network from BA model
  net        <- generate_net(N = 1000 , m = 50 , mode = 1, alpha = 1, s = 0)
  
  net_stats  <- get_statistics(net, only_PA = TRUE)
  result     <- only_A_estimate(net, net_stats)
 
  # plot the estimated attachment function
  plot(result, net_stats)
  
  # true function
  true_A     <- result$estimate_result$center_k
  lines(result$estimate_result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
  
  #### Example 2: a non-log-linear preferential attachment  #########
  # A_k = alpha* log (max(k,1))^beta + 1, with alpha = 2, and beta = 2
  set.seed(1)
  net        <- generate_net(N = 1000 , m = 50 , mode = 3, alpha = 2, beta = 2, s = 0)
  
  net_stats  <- get_statistics(net,only_PA = TRUE)
  result     <- only_A_estimate(net, net_stats)
 
  # plot the estimated attachment function
  plot(result, net_stats)
  
  # true function
  true_A     <- 2 * log(pmax(result$estimate_result$center_k,1))^2 + 1 # true function
  lines(result$estimate_result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
  
  #############################################################################
  #### Example 3: another non-log-linear preferential attachment kernel ############
  set.seed(1)
  # A_k = min(max(k,1),sat_at)^alpha, with alpha = 1, and sat_at = 200
  # inverse variance of the distribution of node fitnesse = 10
  net        <- generate_net(N = 1000 , m = 50 , mode = 2, alpha = 1, sat_at = 200, s = 0)
  net_stats  <- get_statistics(net, only_PA = TRUE)
  
  result     <- only_A_estimate(net, net_stats)
  
  
  # plot the estimated attachment function
  true_A     <- pmin(pmax(result$estimate_result$center_k,1),200)^1 # true function
  plot(result , net_stats, max_A = max(true_A,result$estimate_result$theta))
  lines(result$estimate_result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
  }
}

\concept{preferential attachment}
\concept{attachment function}

