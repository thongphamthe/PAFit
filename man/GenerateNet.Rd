\name{GenerateNet}
\alias{GenerateNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to generate simulated networks from the General Temporal (GT) model}
\description{
This function generates networks from the GT model, a generative temporal network model that includes many well-known models such as the Barabasi-Albert model or the fitness model as special cases. The number of edges of the new node at each time-step can be specified to be fixed, or followed a Poisson distribution.  In the latter case, the mean of the Poisson distribution is either held fixed or increased with time.
}
\usage{
GenerateNet (
 N=1000            , m = 1                , mode = c(1,2,3)  , alpha = 1, 
 beta = 1          , sat_at = 100         , offset = 1       , rate = 0  , 
 shape = 0         , mode_f = c("gamma")  , num_seed = 2     , prob_m = FALSE   ,   
 meanlog = 0       , sdlog = 1            , scale_pareto = 0.6,
 shape_pareto = 2.5, increase = FALSE     , multiple_node = 1, specific_start = NULL, 
 log = FALSE       , custom_PA = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Integer. Total number of nodes in the network (including the nodes in the seed graph). Default value is \eqn{1000}.
}
  \item{m}{
Integer. The number of edges of the new node. Default value is \eqn{1}.
}
\item{mode}{
Integer. Indicates the attachment function to be used in generating the network. If mode == 1, the attachment function is \eqn{A_k = k^\alpha}. If mode == 2, the attachment function is \eqn{A_k = min(k,sat_at)^\alpha}. If mode == 3, the attachment function is \eqn{A_k = \alpha log (k)^\beta}. Default value is \eqn{1}.
} 
  \item{alpha}{
Numeric. If mode == 1, this is the attachment exponent in the attachment function \eqn{A_k = k^\alpha}. If mode == 2, this is the attachment exponenet in the attachment function \eqn{A_k = min(k,sat_at)^\alpha}. If mode == 3, this is the alpha in the attachment function \eqn{A_k = \alpha log (k)^\beta} + 1.
}
\item{beta}{
Numeric. This is the beta in the attachment function \eqn{A_k = \alpha log (k)^\beta} + 1.
}
\item{sat_at}{
Integer. This is the saturation position sat_at in the attachment function \eqn{A_k = min(k,sat_at)^\alpha}.
}

\item{offset}{
Numeric. The attachment value of degree \eqn{0}. Default value is \eqn{1}.
}
  \item{rate}{
Numeric. The rate parameter in the Gamma prior for node fitness. If either rate or shape is \eqn{0}, node fitness \eqn{f} is fixed at \eqn{1} (i.e. Barabasi-Albert model)
}
  \item{shape}{
Numeric. The shape parameter in the Gamma prior for node fitness. If either rate or shape is 0, node fitness \eqn{f} is fixed at \eqn{1} (i.e. Barabasi-Albert model)
}
  \item{mode_f}{
String. Indicates the true distribution for node fitness. "gamma" = gamma distribution, "log_normal" = log-normal distribution. "power_law" = power-law (pareto) distribution. Default value is "gamma".
}
  \item{meanlog}{
Numeric. Mean of the log-normal distribution in log scale. Default value is \eqn{0}.
}
\item{sdlog}{
Numeric. Standard deviation of the log-normal distribution in log scale. Default value is \eqn{1}.
}
\item{scale_pareto}{
Numeric. The scale parameter of the Pareto distribution. Default value is \eqn{0.6}.
}
\item{shape_pareto}{
Numeric. The shape parameter of the Pareto distribution. Default value is \eqn{2.5}.
}
\item{num_seed}{
Integer. The number of nodes of the seed graph (the initial state of the network). Default value is \eqn{2}.
}
 \item{prob_m}{
Logical. Indicates whether we fix the number of edges of a new node as a constant, or let it follows a Poisson 
distribution. If prob_m == TRUE, the number of edges of a new node follows a Poisson distribution. The mean of this distribution depends on the value of increase and log. Default value is FALSE.
}
\item{increase}{
Logical. Indicates whether we increase the mean of the Poisson distribution over time. If increase == FALSE, the mean is fixed at m. If increase == TRUE, the way the mean increases depends on the value of \eqn{log}. Default value is FALSE.
}
\item{multiple_node}{
Numeric. The number of new nodes at each time-step. Default value is \eqn{1}.
}
\item{specific_start}{
Numeric. If \eqn{specific_start} is specified, then all the time-steps from time-step \eqn{1} to \eqn{specific_start} are grouped to become the initial time-step in the final output. This option is used when we want to create a network with a large initial network. Default value is \eqn{NULL}.
}
\item{log}{
Logical. Indicates how to increase the mean of the Poisson distribution. If log == TRUE, the mean increases logarithmically with the number of current nodes. If log == FALSE, the mean increases linearly with the number of current nodes. Default value is FALSE.
}
\item{custom_PA}{
Numeric vector. This is the user-input PA function. If \eqn{custom_PA} is specified, then \eqn{mode} is ignored, and the PA function \eqn{custom_PA} is used to grow the network. Default value is \eqn{NULL}.
}
}

\value{
The output is a List contains the following two fields:
\item{graph}{a 3-column matrix, where each row contains information of one edge, in the form of (from_id, to_id, time_stamp). From_id is the id of the source, to_id is the id of the destination.}
\item{fitness}{fitness values of nodes in the network. The name of each value is the ID of the node.}
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
#Generate a network from the original BA model with alpha = 1, N = 100, m = 1
data   <- GenerateNet(N = 100,m = 1,mode = 1, alpha = 1, shape = 0, rate = 0)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{scale free}
\keyword{BA network}
\keyword{preferential attachment}
\keyword{Bianconi-Barabasi model}
\keyword{ fitness model }% __ONLY ONE__ keyword per line
