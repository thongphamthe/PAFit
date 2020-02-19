\name{joint_estimate}
\alias{joint_estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Joint inference of attachment function and node fitnesses   
}
\description{
  This function jointly estimates the attachment function \eqn{A_k} and node fitnesses \eqn{\eta_i}. It first performs a cross-validation to select the optimal parameters \eqn{r} and \eqn{s}, then estimates \eqn{A_k} and \eqn{eta_i} using that optimal pair with the full data (Ref. 2).
}
\usage{
joint_estimate(net_object                               , 
              net_stat      = get_statistics(net_object), 
              p             = 0.75                      ,
              stop_cond     = 10^-8                     ,
              mode_reg_A    = 0                         , 
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
\item \code{0}: This is the regularization term used in Ref. 1 and 2.  Please refer to Eq. (4) in the tutorial for the definition of the term. It approximately enforces the power-law form \eqn{A_k = k^\alpha}. This is the default value. 
\item \code{1}: Unlike the default, this regularization term exactly enforces the functional form \eqn{A_k = k^\alpha}. Please refer to Eq. (6) in the tutorial for the definition of the term. Its main drawback is it is significantly slower to converge, while its gain over the default one is marginal in most cases.  
}
}
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}

\value{
  Outputs a \code{Full_PAFit_result} object, which is a list containing the following fields:
  \itemize{
    \item \code{cv_data}: a \code{CV_Data} object which contains the cross-validation data. This is the testing data.
    
    \item \code{cv_result}: a \code{CV_Result} object which contains the cross-validation result. Normally the user does not need to pay attention to this data.
    
    \item \code{estimate_result}: this is a \code{PAFit_result} object which contains the estimated attachment function \eqn{A_k}, the estimated fitnesses \eqn{\eta_i} and their confidence intervals. In particular, the important fields are:      
    \itemize{
    \item \code{ratio}: this is the selected value for the hyper-parameter \eqn{r}.
    \item \code{shape}: this is the selected value for the hyper-parameter \eqn{s}.
    \item \code{k} and \code{A}: a degree vector and the estimated PA function.
    \item \code{var_A}: the estimated variance of \eqn{A}.
    \item \code{var_logA}: the estimated variance of \eqn{log A}.
    \item \code{upper_A}: the upper value of the interval of two standard deviations around \eqn{A}.
    \item \code{lower_A}: the lower value of the interval of two standard deviations around \eqn{A}.
    
    \item \code{center_k} and \code{theta}: when we perform binning, these are the centers of the bins and the estimated PA values for those bins. \code{theta} is similar to \code{A} but with duplicated values removed.
     \item \code{var_bin}: the variance of \code{theta}. Same as \code{var_A} but with duplicated values removed.
    \item \code{upper_bin}: the upper value of the interval of two standard deviations around \code{theta}. Same as \code{upper_A} but with duplicated values removed.
    \item \code{lower_bin}: the lower value of the interval of two standard deviations around \code{theta}. Same as \code{lower_A} but with duplicated values removed.
    \item \code{g}: the number of bins used.
    \item \code{alpha} and \code{ci}: \code{alpha} is the estimated attachment exponenet \eqn{\alpha} (when assume \eqn{A_k = k^\alpha}), while \code{ci} is the confidence interval.
    \item \code{loglinear_fit}: this is the fitting result when we estimate \eqn{\alpha}. 
    \item \code{f}: the estimated node fitnesses.
    \item \code{var_f}: the estimated variance of \eqn{\eta_i}.
    \item \code{upper_f}: the estimated upper value of the interval of two standard deviations around \eqn{\eta_i}.
    \item \code{lower_f}: the estimated lower value of the interval of two standard deviations around \eqn{\eta_i}.
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
  
 3. Pham, T., Sheridan, P. & Shimodaira, H. (2020). PAFit: An R Package for the Non-Parametric Estimation of Preferential Attachment and Node Fitness in Temporal Complex Networks. Journal of Statistical Software 92 (3), doi:10.18637/jss.v092.i03. (\url{http://dx.doi.org/10.18637/jss.v092.i03})
}
\seealso{
  See \code{\link{get_statistics}} for how to create summerized statistics needed in this function.
  
  See \code{\link{Jeong}}, \code{\link{Newman}} and \code{\link{only_A_estimate}} for functions to estimate the attachment function in isolation.
  
    See \code{\link{only_F_estimate}} for a function to estimate node fitnesses in isolation.
}

\examples{
\dontrun{
  
  library("PAFit")
  #### Example 1: a linear preferential attachment kernel, i.e., A_k = k ############
  set.seed(1)
  # size of initial network = 100
  # number of new nodes at each time-step = 100
  # Ak = k; inverse variance of the distribution of node fitnesse = 5
  net        <- generate_BB(N        = 1000 , m             = 50 , 
                            num_seed = 100  , multiple_node = 100,
                            s        = 5)
  net_stats  <- get_statistics(net)
  
  # Joint estimation of attachment function Ak and node fitness
  result     <- joint_estimate(net, net_stats)
  
  summary(result)
  
  # plot the estimated attachment function
  true_A     <- pmax(result$estimate_result$center_k,1) # true function
  plot(result , net_stats, max_A = max(true_A,result$estimate_result$theta))
  lines(result$estimate_result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
  
  # plot the estimated node fitnesses and true node fitnesses
  plot(result, net_stats, true = net$fitness, plot = "true_f")
  
  #############################################################################
  #### Example 2: a non-log-linear preferential attachment kernel ############
  set.seed(1)
  # size of initial network = 100
  # number of new nodes at each time-step = 100
  # A_k = alpha* log (max(k,1))^beta + 1, with alpha = 2, and beta = 2
  # inverse variance of the distribution of node fitnesse = 10
  net        <- generate_net(N       = 1000 , m             = 50 , 
                            num_seed = 100  , multiple_node = 100,
                            s        = 10   , mode = 3, alpha = 2, beta = 2)
  net_stats  <- get_statistics(net)
  
  # Joint estimation of attachment function Ak and node fitness
  result     <- joint_estimate(net, net_stats)
  
  summary(result)
  
  # plot the estimated attachment function
  true_A     <- 2 * log(pmax(result$estimate_result$center_k,1))^2 + 1 # true function
  plot(result , net_stats, max_A = max(true_A,result$estimate_result$theta))
  lines(result$estimate_result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
  
  # plot the estimated node fitnesses and true node fitnesses
  plot(result, net_stats, true = net$fitness, plot = "true_f")
  #############################################################################
  #### Example 3: another non-log-linear preferential attachment kernel ############
  set.seed(1)
  # size of initial network = 100
  # number of new nodes at each time-step = 100
  # A_k = min(max(k,1),sat_at)^alpha, with alpha = 1, and sat_at = 100
  # inverse variance of the distribution of node fitnesse = 10
  net        <- generate_net(N       = 1000 , m             = 50 , 
                            num_seed = 100  , multiple_node = 100,
                            s        = 10   , mode = 2, alpha = 1, sat_at = 100)
  net_stats  <- get_statistics(net)
  
  # Joint estimation of attachment function Ak and node fitness
  result     <- joint_estimate(net, net_stats)
  
  summary(result)
  
  # plot the estimated attachment function
  true_A     <- pmin(pmax(result$estimate_result$center_k,1),100)^1 # true function
  plot(result , net_stats, max_A = max(true_A,result$estimate_result$theta))
  lines(result$estimate_result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
  
  # plot the estimated node fitnesses and true node fitnesses
  plot(result, net_stats, true = net$fitness, plot = "true_f")
  }
}

\concept{preferential attachment}
\concept{attachment function}
\concept{fitness}

