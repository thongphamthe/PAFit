\name{generate_net}
\alias{generate_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating networks from preferential attachment and fitness mechanisms}
\description{
This function generates networks from the General Temporal model, a generative temporal network model that includes many well-known models such as the \enc{Erdős–Rényi}{Erdos-Renyi} model, the \enc{Barabási}{Barabasi}-Albert model or the Bianconi-\enc{Barabási}{Barabasi} model as special cases. This function also includes some flexible mechanisms to vary the number of new nodes and new edges at each time-step in order to generate realistic networks. 
}
\usage{
generate_net (N                 = 1000   , 
             num_seed           = 2      , 
             multiple_node      = 1      , 
             specific_start     = NULL   ,
             m                  = 1      ,
             prob_m             = FALSE  ,
             increase           = FALSE  , 
             log                = FALSE  , 
             no_new_node_step   = 0      ,
             m_no_new_node_step = m      ,
             custom_PA          = NULL   ,
             mode               = 1      , 
             alpha              = 1      , 
             beta               = 2      , 
             sat_at             = 100    ,
             offset             = 1      ,
             mode_f             = "gamma", 
             s                  = 10       )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
The parameters can be divided into four groups. 

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
\item{specific_start}{
Positive Integer. If \code{specific_start} is specified, then all the time-steps from time-step \code{1} to \code{specific_start} are grouped to become the initial time-step in the final output. This option is usefull when we want to create a network with a large initial network that follows a scale-free degree distribution. Default value is \code{NULL}.
}
The second group specifies the number of new edges at each time-step:
  \item{m}{
Positive integer. The number of edges of each new node. Default value is \code{1}.
}
 \item{prob_m}{
Logical. Indicates whether we fix the number of edges of each new node as a constant, or let it follows a Poisson distribution. If \code{prob_m == TRUE}, the number of edges of each new node follows a Poisson distribution. The mean of this distribution depends on the value of \code{increase} and \code{log}. Default value is \code{FALSE}.
}
\item{increase}{
Logical. Indicates whether we increase the mean of the Poisson distribution over time. If \code{increase == FALSE}, the mean is fixed at \code{m}. If \code{increase == TRUE}, the way the mean increases depends on the value of \code{log}. Default value is \code{FALSE}.
}
\item{log}{
Logical. Indicates how to increase the mean of the Poisson distribution. If \code{log == TRUE}, the mean increases logarithmically with the number of current nodes. If \code{log == FALSE}, the mean increases linearly with the number of current nodes. Default value is \code{FALSE}.
}
\item{no_new_node_step}{
Non-negative integer. The number of time-steps in which no new node is added, while new edges are added between existing nodes. Default value is \code{0}, i.e., new nodes are always added at each time-step.
}
\item{m_no_new_node_step}{
Positive integer. The number of new edges in the no-new-node steps. Default value is equal to \code{m}. Note that the number of new edges in the no-new-node steps is not effected by the parameters \code{increase} or \code{prob_m}; this number is always the constant specified by \code{m_no_new_node_step}.
}


The third group of parameters specifies the preferential attachment function:
\item{custom_PA}{
Numeric vector. This is the user-input PA function: \eqn{A_0, A_1,..., A_K}. If \code{custom_PA} is specified, then \code{mode} is ignored, and we grow the network using the PA function \code{custom_PA}. Degrees greater than \eqn{K} will have attachment value \eqn{A_k}. Default value is \code{NULL}.
}
\item{mode}{
Integer. Indicates the parametric attachment function to be used in generating the network. If \code{mode == 1}, the attachment function is \eqn{A_k = k^\alpha}. If \code{mode == 2}, the attachment function is \eqn{A_k = min(k,sat.at)^\alpha}. If \code{mode == 3}, the attachment function is \eqn{A_k = \alpha log (k)^\beta}. Default value is \code{1}.
} 
  \item{alpha}{
Numeric. If \code{mode == 1}, this is the attachment exponent in the attachment function \eqn{A_k = k^\alpha}. If \code{mode == 2}, this is the attachment exponenet in the attachment function \eqn{A_k = min(k,sat.at)^\alpha}. If \code{mode == 3}, this is the \eqn{\alpha} in the attachment function \eqn{A_k = \alpha log (k)^\beta + 1}.
}
\item{beta}{
Numeric. This is the beta in the attachment function \eqn{A_k = \alpha log (k)^\beta + 1}.
}
\item{sat_at}{
Integer. This is the saturation position \eqn{sat.at} in the attachment function \eqn{A_k = min(k,sat.at)^\alpha}.
}

\item{offset}{
Numeric. The attachment value of degree \code{0}. Default value is \code{1}.
}
The final group of parameters specifies the distribution from which node fitnesses are generated:
  \item{mode_f}{
String. Possible values:\code{"gamma"}, \code{"log_normal"} or \code{"power_law"}. This parameter indicates the true distribution for node fitness. \code{"gamma"} = gamma distribution, \code{"log_normal"} = log-normal distribution. \code{"power_law"} = power-law (pareto) distribution. Default value is \code{"gamma"}.
}
\item{s}{
Non-negative numeric. The inverse variance parameter. The mean of the distribution is kept at \eqn{1} and the variance is \eqn{1/s} (since node fitnesses are only meaningful up to scale). This is achieved by setting shape and rate parameters of the Gamma distribution to \eqn{s}; setting mean and standard deviation in log-scale of the log-normal distribution to \eqn{-1/2*log (1/s + 1)} and \eqn{(log (1/s + 1))^{0.5}}; and setting shape and scale parameters of the pareto distribution to \eqn{(s+1)^{0.5} + 1} and \eqn{(s+1)^{0.5}/((s+1)^{0.5} + 1)}. If \code{s} is \code{0}, all node fitnesses \eqn{\eta} are fixed at \code{1} (i.e., \enc{Barabási}{Barabasi}-Albert model). The default value is \code{10}.
}
}

\value{
 The output is a \code{PAFit_net} object, which is a List contains the following four fields:
    \item{graph}{a three-column matrix, where each row contains information of one edge, in the form of \code{(from_id, to_id, time_stamp)}. \code{from_id} is the id of the source, \code{to_id} is the id of the destination.}
    \item{type}{a string indicates whether the network is \code{"directed"} or \code{"undirected"}.}
    \item{PA}{a numeric vector contains the true PA function.}
    \item{fitness}{fitness values of nodes in the network. The name of each value is the ID of the node.}
}
\author{
Thong Pham \email{thongpham@thongpham.net}
}
\seealso{
For subsequent estimation procedures, see \code{\link{get_statistics}}.

For simpler functions to generate networks from well-known models, see \code{\link{generate_BA}}, \code{\link{generate_ER}}, \code{\link{generate_BB}} and \code{\link{generate_fit_only}}. }

\examples{
library("PAFit")
#Generate a network from the original BA model with alpha = 1, N = 100, m = 1
net <- generate_net(N = 100,m = 1,mode = 1, alpha = 1, s = 0)
str(net)
plot(net)
}

\concept{fitness model}
\concept{Barabasi-Albert model}
\concept{Bianconi-Barabasi model}
