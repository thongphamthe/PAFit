\name{only_A_estimate}
\alias{only_A_estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimating the attachment function in isolation by PAFit method  
}
\description{
This function estimates the attachment function \eqn{A_k} by PAFit method. It first performs a cross-validation step to select the optimal parameter \eqn{r} for the regularization of \eqn{A_k}, then uses that \eqn{r} to estimate the attachment function. 
}
\usage{
only_A_estimate(raw_net                , 
               net_stat               , 
               stop_cond  = 10^-9     , 
               mode_reg_A = 0         ,
               ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw_net}{
    a three-column matrix that contains the network.
  }
  \item{net_stat}{
    An object of class \code{PAFit_data} which contains summerized statistics needed in estimation. This object is created by the function \code{\link{get_statistics}}.
  }
\item{stop_cond}{Numeric. The iterative algorithm stops when \eqn{abs(h(ii) - h(ii + 1)) / (abs(h(ii)) + 1) < stop.cond} where \eqn{h(ii)} is the value of the objective function at iteration \eqn{ii}. We recommend to choose \code{stop.cond} at most equal to \eqn{10^(- number of digits of h - 2)}, in order to ensure that when the algorithm stops, the increase in posterior probability is less than 1\% of the current posterior probability. Default is \code{10^-9}. This threshold is more than good enough for most applications.}

\item{mode_reg_A}{Binary. Indicates which regularization term is used for \eqn{A_k}:
\itemize{
\item \code{0}: This is the regularization term used in Ref. 1 and 2. It approximately enforces the power-law form \eqn{A_k = k^\alpha}. This is the default value. 
\item \code{1}: Unlike the default, this regularization term exactly enforces the functional form \eqn{A_k = k^\alpha}. Its main drawback is it is significantly slower to converge, while its gain over the default one is marginal in most cases.  
}
}
\item{...}{Other parameters to pass to the internal estimation algorithm.}

}

\value{
  Outputs a list, which contains the following fields.
  \itemize{
  \item \code{cv_data}: a \code{CV_Data} object which contains the cross-validation data. Normally the user does not need to pay attention to this data.
  
 \item \code{cv_result}: a \code{CV_Result} object which contains the cross-validation result. Normally the user does not need to pay attention to this data.
 
 \item \code{estimate_result}: this is a \code{PAFit_result} object which contains the estimated PA function and its confidence interval. It also includes the estimated attachment exponenent \eqn{\alpha} (assuming the model \eqn{A_k = k^\alpha}) in the field \code{alpha}, and the confidence interval of \eqn{\alpha} (in the field \code{ci}) when possible.
}
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).
  
  2. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{www.nature.com/articles/srep32558}).
}
\seealso{
  See \code{\link{get_statistics}} for how to create summerized statistics needed in this function.
  
  See \code{\link{Newman}} and \code{\link{Jeong}} for other methods to estimate the attachment function \eqn{A_k} in isolation.
  
}

\examples{
\dontrun{
  library("PAFit")
  # a network from BA model
  net        <- generate_BB(N = 1000 , m = 1 , mode = 1)
  
  net_stats  <- get_statistics(net$graph)
  result     <- only_A_estimate(net$graph, net_stats)
 
  # plot the estimated attachment function
  plot(result, net_stats)
  
  # true function
  true_A     <- result$estimate_result$center_k
  lines(result$estimate_result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
  }
}

\concept{preferential attachment}
\concept{attachment function}

