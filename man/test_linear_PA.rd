\name{test_linear_PA}
\alias{test_linear_PA}
\title{
 Fitting various distributions to a degree vector
}
\description{
This function implements the method in Handcock and Jones (2004) to fit various distributions to a degree vector. The implemented distributions are Yule, Waring, Poisson, geometric and negative binomial. The Yule and Waring distributions correspond to a preferential attachment situation. In particular, the two distributions correspond to the case of \eqn{A_k = k} for \eqn{k \ge 1} and \eqn{\eta_i = 1} for all \eqn{i} (note that, the number of new edges and new nodes at each time-step are implicitly assumed to be \eqn{1}). 

Thus, if the best fitted distribution, which is chosen by either the Akaike Information Criterion (AIC) or the Bayesian Information Criterion (BIC), is NOT Yule or Waring, then the case of \eqn{A_k = k} for \eqn{k \ge 1} and  \eqn{\eta_i = 1} for all \eqn{i}  is NOT consistent with the observed degree vector.

The method allows the low-tail probabilities to NOT follow the parametric distribution, i.e., \eqn{P(K = k) = \pi_k} for all \eqn{k \le k_min} and \eqn{P(K = k) = f(k,\theta)} for all \eqn{k > k_min}. Here \eqn{k_min} is the degree threshold above which the parametric distribution holds, \eqn{\pi_k} are probabilities of the low-tail, \eqn{f(.,\theta)} is the parametric distribution with parameter vector \eqn{\theta}. 

For fixed \eqn{k_min} and \eqn{f}, \eqn{\pi_k} and \eqn{\theta} can be estimated by Maximum Likelihood Estimation. We can choose the best \eqn{k_min} for each \eqn{f} by comparing the AIC (or BIC). More details can be founded in Handcock and Jones (2004).
}
\usage{
 test_linear_PA(degree_vector)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{degree_vector}{
    a degree vector
  }
}
\value{
  Outputs a \code{Linear_PA_test_result} object which contains the fitting of five distributions to the degree vector: Yule (\code{yule}), Waring (\code{waring}), Poisson (\code{pois}), geometric (\code{geom}) and negative binomial (\code{nb}). In particular, for each distribution, the AIC and BIC are calcualted for each \eqn{k_min}. 
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Handcock MS, Jones JH (2004). “Likelihood-based inference for stochastic models of sexual network formation.” Theoretical Population Biology, 65(4), 413 – 422. ISSN 0040-5809. \url{https://doi.org/10.1016/j.tpb.2003.09.006}. Demography in the 21st Century, \url{http://www.sciencedirect.com/science/article/pii/S0040580904000310}.
}
\examples{
\dontrun{
  library("PAFit")
  set.seed(1)
  net   <- generate_BA(n = 1000)
  stats <- get_statistics(net, only_PA = TRUE)
  u     <- test_linear_PA(stats$final_deg)
  print(u)
}
}

\concept{linear preferential attachment}
\concept{fitting degree distributions}
