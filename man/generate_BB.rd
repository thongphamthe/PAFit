\name{generate_BB}
\alias{generate_BB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulating networks from the Bianconi-\enc{Barabási}{Barabasi} model}
\description{
  This function generates networks from the Bianconi-\enc{Barabási}{Barabasi} model. It is a \sQuote{preferential attachment with fitness} model. In this model, the preferential attachment function is linear, i.e. \eqn{A_k = k}, and node fitnesses are sampled from some probability distribution.  
}
\usage{
generate_BB(N, 
            num_seed       = 2      , 
            multiple_node  = 1      , 
            m              = 1      ,
            mode_f         = "gamma", 
            s              = 10     , 
            meanlog        = 0      , 
            sdlog          = 1      ,
            scale_pareto   = 2      ,
            shape_pareto   = 2       )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  The parameters can be divided into two groups. 
  
The first group specifies basic properties of the network:
    \item{N}{
      Integer. Total number of nodes in the network (including the nodes in the seed graph). Default value is \code{1000}.
    }
  \item{num_seed}{
    Integer. The number of nodes of the seed graph (the initial state of the network). The seed graph is a cycle. Default value is \code{2}.
  }
  \item{multiple_node}{
    Positive integer. The number of new nodes at each time-step. Default value is \code{1}.
  }
\item{m}{
  Positive integer. The number of edges of each new node. Default value is \code{1}.
}

The final group of parameters specifies the distribution from which node fitnesses are generated:
    \item{mode_f}{
      String. Possible values:\code{"gamma"}, \code{"log_normal"} or \code{"power_law"}. This parameter indicates the true distribution for node fitness. \code{"gamma"} = gamma distribution, \code{"log_normal"} = log-normal distribution. \code{"power_law"} = power-law (pareto) distribution. Default value is "gamma".
    }
\item{s}{
    Positive numeric. The inverse variance of the Gamma prior for node fitness. If \code{s} is \code{0}, all node fitnesses \eqn{\eta} are fixed at \code{1} (i.e. \enc{Barabási}{Barabasi}-Albert model)
}
  \item{meanlog}{
    Numeric. Mean of the log-normal distribution in log scale. Default value is \code{0}.
  }
  \item{sdlog}{
    Positive numeric. Standard deviation of the log-normal distribution in log scale. Default value is \code{1}.
  }
  \item{scale_pareto}{
    Numeric. The scale parameter of the Pareto distribution. Default value is \code{2}.
  }
  \item{shape_pareto}{
    Numeric. The shape parameter of the Pareto distribution. Default value is \code{2}.
  }
}

\value{
  The output is a List contains the following two fields:
    \item{graph}{a three-column matrix, where each row contains information of one edge, in the form of \code{(from_id, to_id, time_stamp)}. \code{from_id} is the id of the source, \code{to_id} is the id of the destination.}
  \item{fitness}{fitness values of nodes in the network. The name of each value is the ID of the node.}
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Bianconni, G. & \enc{Barabási}{Barabasi}, A. (2001). Competition and multiscaling in evolving networks. Europhys. Lett., 54, 436 (\url{http://iopscience.iop.org/article/10.1209/epl/i2001-00260-6/meta}).
}
\seealso{
  For subsequent estimation procedures, see \code{\link{get_statistics}}.
  
  For other functions to generate networks, see \code{\link{generate_net}}, \code{\link{generate_BA}}, \code{\link{generate_ER}} and \code{\link{generate_fit_only}}. }

\examples{
  library("PAFit")
  # Generate a network from the BB model with alpha = 1, N = 100, m = 1
  # The inverse variance of distribution of node fitnesses is s = 10
  net <- generate_BB(N = 100,m = 1,mode = 1, s = 10)
  str(net)
}
\concept{Bianconi-\enc{Barabási}{Barabasi} model}
