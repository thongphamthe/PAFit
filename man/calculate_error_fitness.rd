\name{calculate_error_fitness}
\alias{calculate_error_fitness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A function to calculate the relative error between the true and estimated node fitnesses.
}
\description{
  This function calculates the relative error \eqn{e_\eta} between the true node fitnesses \eqn{\eta} and the estimated node fitnesses \eqn{\hat{eta}}. \eqn{e_\eta} is defined as \eqn{e_\eta = 1/{N}\sum_{i = 1}^{N} (\eta_i - \hat{\eta}_i)^2/\eta_i^2}. 
}
\usage{
  calculate_error_fitness(true, estimate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{true}{
    Vector. The true node fitnesses.
  }
  \item{estimate}{
    Vector. The estimated node fitnesses. 
  }
}
\value{
  Outputs the relative error \eqn{e_f}.
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
  net        <- GenerateNet(N = 1000 , m = 1 , mode = 1 , alpha = 1 , shape = 100, rate = 100)
  net_stats  <- GetStatistics(net$graph)
  result     <- PAFit(net_stats)
  error_f    <- calculate_error_fitness(true = net$fitness, estimate = result$f)
  print(error_f)
}
