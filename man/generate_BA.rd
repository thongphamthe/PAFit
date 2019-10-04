\name{generate_BA}
\alias{generate_BA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulating networks from the generalized Barabasi-Albert model}
\description{
  This function generates networks from the generalized \enc{Barabási}{Barabasi}-Albert model. In this model, the preferential attachment function is power-law, i.e. \eqn{A_k = k^\alpha}, and node fitnesses are all equal to \eqn{1}. It is a wrapper of the more powerful function \code{\link{generate_net}}. 
}
\usage{
generate_BA(N              = 1000, 
            num_seed       = 2   , 
            multiple_node  = 1   , 
            m              = 1   ,
            alpha          = 1)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
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
\item{alpha}{
Numeric. This is the attachment exponent in the attachment function \eqn{A_k = k^\alpha}. }
}


\value{
   The output is a \code{PAFit_net} object, which is a List contains the following four fields:
    \item{graph}{a three-column matrix, where each row contains information of one edge, in the form of \code{(from_id, to_id, time_stamp)}. \code{from_id} is the id of the source, \code{to_id} is the id of the destination.}
    \item{type}{a string indicates whether the network is \code{"directed"} or \code{"undirected"}.}
    \item{PA}{a numeric vector contains the true PA function.}
  \item{fitness}{fitness values of nodes in the network. The fitnesses are all equal to \eqn{1}.}
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Albert, R. & \enc{Barabási}{Barabasi}, A. (1999). Emergence of scaling in random networks. Science, 286,509–512 (\url{http://science.sciencemag.org/content/286/5439/509}).
}
\seealso{
For subsequent estimation procedures, see \code{\link{get_statistics}}.

For other functions to generate networks, see \code{\link{generate_net}}, \code{\link{generate_ER}}, \code{\link{generate_BB}} and \code{\link{generate_fit_only}}. }

\examples{
  library("PAFit")
  # generate a network from the BA model with alpha = 1, N = 100, m = 1
  net <- generate_BA(N = 100)
  str(net)
  plot(net)
}

\concept{Barabasi-Albert model}
