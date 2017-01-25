\name{PAFit-package}
\alias{PAFit-package}
\alias{FitMultinomial}
\alias{calculate_error_PA}
\alias{calculate_error_fitness}
\docType{package} 
\title{A statistical sound method for joint estimation of the attachment function and node fitness in a temporal complex network by maximizing a suitable penalized log-likelihood function}
\description{
A method for estimating jointly the attachment function \eqn{A_k} and node fitness \eqn{\eta_i} of a temporal complex network is implemented in this package.  The network's growth is assumed to follow a modified version of the fitness model, in which newly added edges, including newly edges between existed nodes, are connected to a degree \eqn{k} node \eqn{v_i} with probability proportional to the product of the attachment value \eqn{A_k} and the fitness value \eqn{\eta_i}.  The method makes no assumption on the functional form of either \eqn{A_k} or \eqn{\eta_i}.  By choosing suitable regularizations, good estimations of \eqn{A_k} and \eqn{\eta_i} can be obtained by maximizing the corresponding penalized log-likelihood function.  We also implement a fast estimation of confidence intervals based on the Hessian of the penalized log likelihood.  See the accompanying vignette for a tutorial.

For a list of references, please run the command: citation("PAFit").
}
\details{
\tabular{ll}{
Package: \tab PAFit\cr
Type: \tab Package\cr
Version: \tab 0.9.3\cr
Date: \tab 2017-01-20\cr
License: \tab GPL-3 \cr
}
\itemize{
\item PAFit: estimates the Preferential Attachment function and fitness function in a temporal complex network.
\item GenerateNet: generates simulated networks based on the Barabasi-Albert model or the fitness model.
\item GetStatistics: summarizes a matrix of edges into summary statistics ready for applying the PAFit function.
}
}
\author{
Thong Pham, Paul Sheridan, Hidetoshi Shimodaira.
Maintainer: Thong Pham \email{thongpham@thongpham.net}
}


\references{
1. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Nonparametric Estimation of the Preferential Attachment Function in Complex Networks: Evidence of Deviations from Log Linearity, Proceedings of ECCS 2014, 141-153 (Springer International Publishing) (\url{http://dx.doi.org/10.1007/978-3-319-29228-1_13}).

2. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).

3. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{www.nature.com/articles/srep32558}).
}


\keyword{ package }
\keyword{fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{Barabasi-Albert model}
\keyword{MM algorithims}
\keyword{Bianconi-Barabasi model}
\keyword{fitness}
\keyword{maximum likelihood function}
\examples{
library("PAFit")
data   <- GenerateNet(N = 100,m = 1,alpha = 1, mode = 1, shape = 0, rate = 0)
stats  <- GetStatistics(data$graph)
result <- PAFit(stats,only_PA = TRUE,stop_cond = 10^-3)
plot(x = result,data = stats,plot = "A")
}
