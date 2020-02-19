\name{PAFit-package}
\alias{PAFit-package}
\alias{PAFit}
\docType{package} 
\title{Generative Mechanism Estimation in Temporal Complex Networks}
\description{
A package for estimating preferential attachment and node fitness generative mechanisms in temporal complex networks. References: Thong Pham et al. (2015) <10.1371/journal.pone.0137796>, Thong Pham et al. (2016) <doi:10.1038/srep32558>, Thong Pham et al. (2020) <doi:10.18637/jss.v092.i03>.}
\details{
\tabular{ll}{
Package: \tab PAFit\cr
Type: \tab Package\cr
Version: \tab 1.0.1.8\cr
Authors: \tab Thong Pham, Paul Sheridan, Hidetoshi Shimodaira\cr
Maintainer: \tab Thong Pham \email{thongphamthe@gmail.com}\cr
Date: \tab 2019-02-17\cr
License: \tab GPL-3 \cr
}
The PAFit package provides a comprehensive framework to deal with growth mechanisms of temporal complex networks. In particular, it implements functions to simulate various temporal network models, gather essential network statistics from raw input data, and use these summarized statistics in the estimation of the attachment function \eqn{A_k} and node fitnesses \eqn{\eta_i}. The heavy computational parts of the package are implemented in \code{C++} through the use of the Rcpp package. Furthermore, users with a multi-core machine can enjoy a hassle-free speed up through OpenMP parallelization mechanisms implemented in the code. Apart from the main functions, the package also includes a real-world collaboration network dataset between scientists in the field of complex networks (\code{\link{coauthor.net}}). The main package functionalities are as follows. 

Firstly, most well-known temporal network models based on the preferential attachment (PA) and node fitness mechanisms can be easily simulated using the package. PAFit implements \code{\link{generate_BA}} for the \enc{Barabási}{Barabasi}-Albert (BA) model, \code{\link{generate_ER}} for the growing \enc{Erdős–Rényi}{Erdos-Renyi} (ER) model, \code{\link{generate_BB}} for the Bianconi-\enc{Barabási}{Barabasi} (BB) model and \code{\link{generate_fit_only}} for the Caldarelli model. These functions have many customizable options, for example the number of new edges at each time-step are tunable stochastic variables. They are actually wrappers of the more powerful \code{\link{generate_net}} function, which simulates networks with more flexible attachment function and node fitness settings.

Secondly, the function \code{\link{get_statistics}} efficiently collects all temporal network summary statistics. We note that \code{\link{get_statistics}} automatically handles both directed and undirected networks. It returns a list containing many statistics that can be used to characterize the network growth process. Notable fields are \code{m_tk} containing the number of new edges that connect to a degree-\eqn{k} node at time-step \eqn{t}, and \code{node_degree} containing the degree sequence, i.e., the degree of each node at each time-step.

The most important functionality of the package is estimating the attachment function and node fitnesses of a temporal network. This is implemented through various methods. There are three usages: estimation of the attachment function in isolation, estimation of the node fitnesses in isolation, and the joint estimation of the attachment function and node fitnesses. 
\itemize{
\item The functions for estimating the attachment function in isolation are: \code{\link{Jeong}} for Jeong's method (Ref. 1), \code{\link{Newman}} for Newman's method (Ref. 2), and \code{\link{only_A_estimate}} for the PAFit method (Ref. 3). 
\item For estimation of node fitnesses in isolation, \code{\link{only_F_estimate}} implements a variant of the PAFit method (Ref. 4). 
\item For the joint estimation of the attachment function and node fitnesses, we implement the full version of the PAFit method in \code{\link{joint_estimate}} (Ref. 4). 
}
In all cases, the input of these functions is the output object of the function \code{\link{get_statistics}}. The output object of these functions contains the estimation results as well as some additional information pertaining to the estimation process. The estimated attachment function and/or node fitnesses can be plotted by using the \code{plot} command directly on this output object. This will visualize not only the estimated results but also the remaining uncertainties when possible.

}

\author{
Thong Pham \email{thongpham@thongpham.net}, Paul Sheridan, and Hidetoshi Shimodaira. 
}

\references{
1. Jeong, H., \enc{Néda}{Neda}, Z. & \enc{Barabási}{Barabasi}, A. (2003). Measuring Preferential Attachment in Evolving Networks. \emph{Europhysics Letters} 61(61):567-572. doi:10.1209/epl/i2003-00166-9 (\url{https://iopscience.iop.org/article/10.1209/epl/i2003-00166-9}).

2. Newman, M. (2001). Clustering and Preferential Attachment in Growing Networks. \emph{Physical Review E} 64(2):025102. doi:10.1103/PhysRevE.64.025102 (\url{https://doi.org/10.1209/epl/i2003-00166-9}).

3. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. \emph{PLOS ONE} 10(9):e0137796. doi:10.1371/journal.pone.0137796 (\url{https://dx.doi.org/10.1371/journal.pone.0137796}).

4. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. \emph{Scientific Reports} 6, Article number: 32558. doi:10.1038/srep32558 (\url{http://www.nature.com/articles/srep32558}).

5. Pham, T., Sheridan, P. & Shimodaira, H. (2020). PAFit: An R Package for the Non-Parametric Estimation of Preferential Attachment and Node Fitness in Temporal Complex Networks. \emph{Journal of Statistical Software} 92 (3), doi:10.18637/jss.v092.i03. (\url{http://dx.doi.org/10.18637/jss.v092.i03})
}

\seealso{

 See the acompanying vignette for a tutorial.

 See also the \href{https://github.com/thongphamthe/PAFit}{GitHub page}.

}


\concept{temporal complex networks}
\concept{preferential attachment}
\concept{attachment function}
\concept{fitness model}
\concept{Jeong's method}
\concept{Corrected Newman's method}
\concept{Barabasi-Albert model}
\concept{Bianconi-Barabasi model}

\examples{
\dontrun{
  ### Jointly estimate the attachment function and node fitnesses
   library("PAFit")
   set.seed(1)
  # a Bianconi-Barabasi network 
  # size of initial network = 100
  # number of new nodes at each time-step = 100
  # Ak = k; inverse variance of distribution of fitness: s = 10
  net        <- generate_BB(N        = 1000 , m             = 10 , 
                            num_seed = 100  , multiple_node = 100,
                            s        = 10)
  net_stats  <- get_statistics(net)
  
  #Joint estimation of attachment function Ak and node fitness
  result     <- joint_estimate(net, net_stats)
  
  summary(result)
  
  # plot the estimated attachment function
  plot(result, net_stats)
  
  # true function
  true_A     <- pmax(result$estimate_result$center_k,1)
  lines(result$estimate_result$center_k, true_A, col = "red") # true line
  legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")
  #plot distribution of estimated node fitnesses
  plot(result, net_stats, plot = "f")
  
  #plot the estimated node fitnesses and true node fitnesses
  plot(result, net_stats, true = net$fitness, plot = "true_f")
}
}
