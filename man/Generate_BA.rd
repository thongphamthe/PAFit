\name{Generate_BA}
\alias{Generate_BA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulating networks from the generalized \enc{Barabási}{Barabasi}-Albert model}
\description{
  This function generates networks from the generalized \enc{Barabási}{Barabasi}-Albert model. In this model, the preferential attachment function is power-law, i.e. \eqn{A_k = k^\alpha}, and node fitnesses are all equal to \eqn{1}. It is a wrapper of the more powerful function \code{\link{GenerateNet}}. 
}
\usage{
Generate_BA(N                  , 
            num_seed       = 2 , 
            multiple_node  = 1 , 
            m              = 1 ,
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
  The output is a List contains the following two fields:
    \item{graph}{a three-column matrix, where each row contains information of one edge, in the form of \code{(from_id, to_id, time_stamp)}. \code{from_id} is the id of the source, \code{to_id} is the id of the destination.}
  \item{fitness}{fitness values of nodes in the network. The fitnesses are all equal to \eqn{1}.}
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}
\references{
  1. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Nonparametric Estimation of the Preferential Attachment Function in Complex Networks: Evidence of Deviations from Log Linearity, Proceedings of ECCS 2014, 141-153 (Springer International Publishing) (\url{http://dx.doi.org/10.1007/978-3-319-29228-1_13}).
  
  2. Pham, T., Sheridan, P. & Shimodaira, H. (2015). PAFit: A Statistical Method for Measuring Preferential Attachment in Temporal Complex Networks. PLoS ONE 10(9): e0137796. doi:10.1371/journal.pone.0137796 (\url{http://dx.doi.org/10.1371/journal.pone.0137796}).
  
  3. Pham, T., Sheridan, P. & Shimodaira, H. (2016). Joint Estimation of Preferential Attachment and Node Fitness in Growing Complex Networks. Scientific Reports 6, Article number: 32558. doi:10.1038/srep32558   (\url{www.nature.com/articles/srep32558}).
}
\seealso{
For subsequent estimation procedures, see \code{\link{GetStatistics}} and \code{\link{CreateDataCV}}.

For other functions to generate networks, see \code{\link{GenerateNet}}, \code{\link{Generate_ER}}, \code{\link{Generate_BB}} and \code{\link{Generate_fitonly}}. }

\examples{
  library("PAFit")
  #Generate a network from the original BA model with alpha = 1, N = 100, m = 1
  net <- Generate_BA(N = 100)
  str(net)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{scale free}
\keyword{BA network}
\keyword{preferential attachment}
\keyword{Bianconi-Barabasi model}
\keyword{ fitness model }% __ONLY ONE__ keyword per line
