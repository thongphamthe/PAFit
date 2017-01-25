\name{GetStatistics}
\alias{GetStatistics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to summarize input data into sufficient statistics for estimating the attachment function and node fitness
}
\description{
The function summarizes input data into sufficient statistics for estimating the attachment function and node fitness, together with additional information about the data, such as the total number of nodes, the number of time-steps, the maximum degree, and the final degree of the network, etc. .  It also provides mechanisms to automatically deal with very large datasets by binning the degree, setting a degree threshold, or discarding time-steps. 
}
\usage{GetStatistics(data,net_type = c("directed","undirected"), only_PA = FALSE, 
                     Binning = TRUE, G = 200, start_deg = 0,  
                     deg_threshold = 0, CompressMode = c(0,1,2,3), 
                     CompressRatio = 0.5 , CustomTime = NULL, 
                     only_true_deg_matrix = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix. A 3-column matrix whose each row contains information of one edge in the form (from_node id, to_node id, time_stamp). from_node id is the id of the source node. to_node id is the id of the destination node. time_stamp is the arrival time of the edge. from_node id and to_node id are assumed to be integers starting from \eqn{0}. time_stamp can be either numeric or string. The value of a time-stamp can be arbitrary, but we assume that a smaller time_stamp (regarded so by the \eqn{sort} function in \eqn{R}) represents an earlier arrival time.}
\item{net_type}{
String. Indicates the type of network. Default value is \eqn{"directed"}.
}
\item{only_PA}{
Logical. Indicates whether the statistics for estimating \eqn{A_k} are summarized (if TRUE, saving memory at the cost of unable to estimate node fitness). Default value is \eqn{FALSE}. 
}
\item{Binning}{
Logical. Indicates whether the degree should be binned together. Default value is \eqn{TRUE}.
}
\item{G}{
Integer. Number of bins. Default value is \eqn{200}.
}
 \item{start_deg}{
Integer. The degree from which the program start to binning the degree together. Default value is \eqn{0}.
}
\item{deg_threshold}{
Integer. We only estimate the fitnesses of nodes whose number of new edges acquired is at least deg_threshold. The fitnesses of all other nodes are fixed at \eqn{1}. Default value is \eqn{0}.
}
\item{CompressMode}{
Integer. Indicates whether the timeline should be compressed. The value of CompressMode:

0: No compression

1: Compressed by using a subset of time-steps. The time stamps in this subset are equally spaced. The size of this subset is \eqn{CompressRatio} times the size of set of all time stamps.

2: Compressed by only starting from the first time-step when \eqn{CompressRatio*100} percentages of the total number of edges (in the final state of the network) had already been added to the network.

3: This mode offers the most flexibility, but requires user to supply the time stamps in \eqn{CustomTime}. Only time stamps in this \eqn{CustomTime} will be used. This mode can be used, for example, when investigating the change of the attachment function or node fitness in different time intervals.

 Default value is \eqn{0}, i.e. no compression.
}
  \item{CompressRatio}{
Numeric. Indicates how much we should compress if CompressMode is 1 or 2. Default value is \eqn{0.5}.
}
  \item{CustomTime}{
Vector. Custom time stamps. Only effective if CompressMode == 3. In that case, only these time stamps are used.
}
\item{only_true_deg_matrix}{
Logical. Return only the true degree matrix (without binning). The motivation for this option is that sometimes we only want to get a degree matrix that summerizes the growth process of the whole network for plotting etc. Default value is \eqn{FALSE}.
}
}
\value{An object of class PAFitData, which is a list. Some important fields are:
\item{offset_tk}{A matrix where the (t,k+1) element is the number of nodes with degree k at time t, counting among all the nodes whose number of new edges acquired is less than \eqn{deg_thresh}}
\item{n_tk}{A matrix where the (t,k+1) element is the number of nodes with degree \eqn{k} at time \eqn{t}}
\item{m_tk}{A matrix where the (t,k+1) element is the number of new edges connect to a degree-\eqn{k} node at time \eqn{t}}
\item{Sum_m_k}{A vector where the (k+1)-th element is the total number of edges that linked to a degree \eqn{k} node, counting over all time steps}
\item{node_degree}{A matrix recording the degree of all nodes (that satisty degree_threshold) at each time step}
\item{m_t}{A vector where the t-th element is the number of new edges at time t}
\item{z_j}{A vector where the j-th element is the total number of edges that linked to node j}
\item{N}{Numeric. The number of nodes in the network}
\item{T}{Numeric. The number of time steps}
\item{deg.max}{Numeric. The maximum degree in the final network}
\item{node_id}{A vector contains the id of all nodes}
\item{final_deg}{A vector contains the final degree of all nodes (including those that do not satisfy degree_threshold)}
\item{deg_thresh}{Numeric. The specified degree threshold.}
\item{f_position}{Numeric vector. The index in the node_id vector of the nodes we want to estimate (i.e. nodes whose number of new edges acquired is at least deg_thresh)}
\item{start_deg}{Numeric. The degree at which we start binning.}
\item{begin_deg}{Numeric vector contains the beginning degree of each bin}
\item{end_deg}{Numeric vector contains the ending degree of each bin}
\item{interval_length}{Numeric vector contains the length of each bin.}
\item{Binning}{Logical. Indicates whether binning was applied or not.}
\item{G}{Integer. Number of bins}
\item{TimeCompressMode}{Integer. The mode of time compression.}
\item{T_compressed}{Integer. The number of time stamps actually used}
\item{compressed_unique_time}{The time stamps that are actually used}
\item{CompressRatio}{Numeric.}
\item{CustomTime}{Vector. The time stamps specified by user.}
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
data   <- GenerateNet(N = 100,m = 1,mode = 1, alpha = 1, shape = 5, rate = 5)
stats  <- GetStatistics(data$graph)
}

\keyword{fitness model}
\keyword{Preferential Attachment function}
\keyword{Minorize-Maximization algorithms}
\keyword{Bianconi-Barabasi model}
\keyword{Barabasi-Albert model}% __ONLY ONE__ keyword per line
