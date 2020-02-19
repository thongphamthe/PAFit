\name{only_F_estimate}
\alias{only_F_estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimating node fitnesses in isolation   
}
\description{
  This function estimates node fitnesses \eqn{\eta_i} assusming either \eqn{A_k = k} (i.e. linear preferential attachment) or \eqn{A_k = 1} (i.e. no preferential attachment). The method has a hyper-parameter \eqn{s}. It first performs a cross-validation to select the optimal parameter \eqn{s} for the prior of \eqn{\eta_i}, then estimates \eqn{eta_i} with the full data (Ref. 1).
}
\usage{
only_F_estimate(net_object                             , 
               net_stat    = get_statistics(net_object), 
               p           = 0.75                      ,
               stop_cond   = 10^-8                     , 
               model_A     = "Linear"                  ,
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
  
  \item{model_A}{String. Indicates which attachment function \eqn{A_k} we assume:
      \itemize{
        \item \code{"Linear"}: We assume \eqn{A_k = k}, i.e. the Bianconi-\enc{Barabási}{Barabasi} model (Ref. 2).
        \item \code{"Constant"}: We assume \eqn{A_k = 1}, i.e. the Caldarelli model (Ref. 3).
      }
  }
    \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
  }
}

\value{
   Outputs a \code{Full_PAFit_result} object, which is a list containing the following fields:
  \itemize{
    \item \code{cv_data}: a \code{CV_Data} object which contains the cross-validation data. Normally the user does not need to pay attention to this data.
    
    \item \code{cv_result}: a \code{CV_Result} object which contains the cross-validation result. Normally the user does not need to pay attention to this data.
    
    \item \code{estimate_result}: this is a \code{PAFit_result} object which contains the estimated node fitnesses and their confidence intervals. In particular, the important fields are:      
    \itemize{
    \item \code{shape}: this is the selected value for the hyper-parameter \eqn{s}.
    \item \code{g}: the number of bins used.
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
  1. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{http://www.nature.com/articles/srep32558}).
  
  2. Bianconni, G. & \enc{Barabási}{Barabasi}, A. (2001). Competition and multiscaling in evolving networks. Europhys. Lett., 54, 436 (\url{http://iopscience.iop.org/article/10.1209/epl/i2001-00260-6/meta}).
  
  3. Caldarelli, G., Capocci, A. , De Los Rios, P. & \enc{Muñoz}{Munoz}, M.A. (2002). Scale-Free Networks from Varying Vertex Intrinsic Fitness. Phys. Rev. Lett., 89, 258702 (\url{http://link.aps.org/doi/10.1103/PhysRevLett.89.258702}).
  
}
\seealso{
  See \code{\link{get_statistics}} for how to create summerized statistics needed in this function.
  
  See \code{\link{joint_estimate}} for the method to jointly estimate the attachment function \eqn{A_k} and node fitnesses \eqn{\eta_i}.
  
}

\examples{
\dontrun{
  library("PAFit")
  set.seed(1)
  # size of initial network = 100
  # number of new nodes at each time-step = 100
  # Ak = k; inverse variance of the distribution of node fitnesse = 10
  net        <- generate_BB(N        = 1000 , m             = 50 , 
                            num_seed = 100  , multiple_node = 100,
                            s        = 10)
                            
  net_stats  <- get_statistics(net)
  
  # estimate node fitnesses in isolation, assuming Ak = k
  result     <- only_F_estimate(net, net_stats)
 
  # plot the estimated node fitnesses and true node fitnesses
  plot(result, net_stats, true = net$fitness, plot = "true_f")
  }
}

\concept{fitness model}