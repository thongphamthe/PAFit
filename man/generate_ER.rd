\name{generate_ER}
\alias{generate_ER}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulating networks from the Erdos-Renyi model}
\description{
  This function generates networks from the \enc{Erdős–Rényi}{Erdos-Renyi} model. In this model, the preferential attachment function is a constant function, i.e. \eqn{A_k = 1}, and node fitnesses are all equal to \eqn{1}. It is a wrapper of the more powerful function \code{\link{generate_net}}.
}
\usage{
  generate_ER(N              = 1000, 
              num_seed       = 2   , 
              multiple_node  = 1   , 
              m              = 1)
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
  1. \enc{Erdös}{Erdos} P. & \enc{Rényi}{Renyi} A.. On random graphs. Publicationes Mathematicae Debrecen. 1959;6:290–297 (\url{http://snap.stanford.edu/class/cs224w-readings/erdos59random.pdf}).
}
\seealso{
For subsequent estimation procedures, see \code{\link{get_statistics}}.

For other functions to generate networks, see \code{\link{generate_net}}, \code{\link{generate_BA}}, \code{\link{generate_BB}} and \code{\link{generate_fit_only}}. }

\examples{
  library("PAFit")
  # generate a network from the ER model with N = 1000 nodes
  net <- generate_ER(N = 1000)
  str(net)
  plot(net)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ER model}
