\name{calculate_error_PA}
\alias{calculate_error_PA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculating Preferential Attachment Estimation Error
}
\description{
  This function calculates the relative error \eqn{e_A} between the true attachment function \eqn{A_k} and the estimated attachment function \eqn{\hat{A}_k}. \eqn{e_A} is defined as \eqn{e_A = 1/{K-start_deg + 1}\sum_{k = start_deg}^{K} (A_k - \hat{A}_k)^2/A^2_k}. 
}
\usage{
calculate_error_PA(k           ,  A        , start_deg = 0 ,
                   mode   = 1  , alpha = 1 , beta      = 1 , 
                   sat_at = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
    Integer vector. The degree sequence at which \eqn{A_k} are estimated.
  }
  \item{A}{
    Vector. The estimated attachment function. 
  }
  \item{start_deg}{Integer. The starting degree from which the relative error \eqn{e_A} is calculated. Default is \eqn{0}.}
 \item{mode}{
Integer. Indicates the attachment function to be used in generating the network. If \code{mode == 1}, the attachment function is \eqn{A_k = k^\alpha}. If \code{mode == 2}, the attachment function is \eqn{A_k = min(k,sat_at)^\alpha}. If \code{mode == 3}, the attachment function is \eqn{A_k = \alpha log (k)^\beta}. Default value is \eqn{1}.
} 
  \item{alpha}{
Numeric. If \code{mode == 1}, this is the attachment exponent in the attachment function \eqn{A_k = k^\alpha}. If \code{mode == 2}, this is the attachment exponenet in the attachment function \eqn{A_k = min(k,sat_at)^\alpha}. If \code{mode == 3}, this is the alpha in the attachment function \eqn{A_k = \alpha log (k)^\beta} + 1.
}
\item{beta}{
Numeric. This is the beta in the attachment function \eqn{A_k = \alpha log (k)^\beta} + 1.
}
\item{sat_at}{
Integer. This is the saturation position \eqn{sat_at} in the attachment function \eqn{A_k = min(k,sat_at)^\alpha}.
}
}
\value{
  Outputs the relative error \eqn{e_A}.
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
  error_A    <- calculate_error_PA(k = result$k , A = result$A , mode = 1 , alpha = 1)
  print(error_A)
}
