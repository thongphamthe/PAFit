\name{get_statistics}
\alias{get_statistics}
\alias{PAFit_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Getting summarized statistics from input data  
}
\description{
The function summarizes input data into sufficient statistics for estimating the attachment function and node fitness, together with additional information about the data, such as total number of nodes, number of time-steps, maximum degree, and the final degree of the network, etc. .  It also provides mechanisms to automatically deal with very large datasets by binning the degree, setting a degree threshold, or grouping time-steps. 
}
\usage{
get_statistics(net_object, only_PA  = FALSE , 
               only_true_deg_matrix = FALSE ,
               binning              = TRUE  , g              = 50    , 
               deg_threshold        = 0     , 
               compress_mode        = 0     , compress_ratio = 0.5   , 
               custom_time          = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

The parameters can be divided into four groups. The first group specifies input data and how the data will be summerized:

  \item{net_object}{
An object of class \code{PAFit_net}. You can use the function \code{\link{as.PAFit_net}} to convert from an edgelist matrix, function \code{\link{from_igraph}} to convert from an \code{igraph} object, function \code{\link{from_networkDynamic}} to convert from a \code{networkDynamic} object, and function \code{\link{graph_from_file}} to read from a file.}

\item{only_PA}{
Logical. Indicates whether only the statistics for estimating \eqn{A_k} are summarized. if \code{TRUE}, the statistics for estimating \eqn{\eta_i} are NOT collected. This will save memory at the cost of unable to estimate node fitness). Default value is \code{FALSE}. 
}
\item{only_true_deg_matrix}{
Logical. Return only the true degree matrix (without binning), and no other statistics is returned. The result cannot be used in \code{PAFit} function to estimate PA or fitness. The motivation for this option is that sometimes we only want to get a degree matrix that summerizes the growth process of a very big network for plotting etc. Default value is \code{FALSE}.
}
Second group of parameters specifies how to bin the degrees:
\item{binning}{
Logical. Indicates whether the degree should be binned together. Default value is \code{TRUE}.
}
\item{g}{
Positive integer. Number of bins. Should be at least \code{3}. Default value is \code{50}.
}
Third group contains a single parameter specifying how to reduce the number of node fitnesses:
\item{deg_threshold}{
Integer. We only estimate the fitnesses of nodes whose number of new edges acquired is at least \code{deg_threshold}. The fitnesses of all other nodes are fixed at \code{1}. Default value is \code{0}.
}
Last group of parameters specifies how to group the time-stamps: 
\item{compress_mode}{
Integer. Indicates whether the timeline should be compressed. The value of CompressMode:

\code{0}: No compression

\code{1}: Compressed by using a subset of time-steps. The time stamps in this subset are equally spaced. The size of this subset is \code{CompressRatio} times the size of the set of all time stamps.

\code{2}: Compressed by only starting from the first time-step when \eqn{CompressRatio*100} percentages of the total number of edges (in the final state of the network) had already been added to the network.

\code{3}: This mode offers the most flexibility, but requires user to supply the time stamps in \code{CustomTime}. Only time stamps in this \code{CustomTime} will be used. This mode can be used, for example, when investigating the change of the attachment function or node fitness in different time intervals.

 Default value is \code{0}, i.e. no compression.
}
  \item{compress_ratio}{
Numeric. Indicates how much we should compress if CompressMode is \code{1} or {2}. Default value is \code{0.5}.
}
  \item{custom_time}{
Vector. Custom time stamps. This vector is a subset of the vector that contains all time-stamps. Only effective if \code{CompressMode == 3}. In that case, only these time stamps are used.
}
}
\value{An object of class \code{PAFit_data}, which is a list. Some important fields are:
\item{offset_tk}{A matrix where the \code{(t,k+1)} element is the number of nodes with degree \eqn{k} at time \eqn{t}, counting among all the nodes whose number of new edges acquired is less than \code{deg_thresh}}
\item{n_tk}{A matrix where the \code{(t,k+1)} element is the number of nodes with degree \eqn{k} at time \eqn{t}}
\item{m_tk}{A matrix where the \code{(t,k+1)} element is the number of new edges connect to a degree-\eqn{k} node at time \eqn{t}}
\item{sum_m_k}{A vector where the \code{(k+1)}-th element is the total number of edges that linked to a degree \eqn{k} node, counting over all time steps}
\item{node_degree}{A matrix recording the degree of all nodes (that satisty \code{degree_threshold} condition) at each time step}
\item{m_t}{A vector where the \code{t}-th element is the number of new edges at time \eqn{t}}
\item{z_j}{A vector where the \code{j}-th element is the total number of edges that linked to node \eqn{j}}
\item{N}{Numeric. The number of nodes in the network}
\item{T}{Numeric. The number of time steps}
\item{deg_max}{Numeric. The maximum degree in the final network}
\item{node_id}{A vector contains the id of all nodes}
\item{final_deg}{A vector contains the final degree of all nodes (including those that do not satisfy the \code{degree_threshold} condition)}
\item{deg_thresh}{Integer. The specified degree threshold.}
\item{f_position}{Numeric vector. The index in the \code{node_id} vector of the nodes we want to estimate (i.e. nodes whose number of new edges acquired is at least \code{deg_thresh})}
\item{start_deg}{Integer. The specified degree at which we start binning.}
\item{begin_deg}{Numeric vector contains the beginning degree of each bin}
\item{end_deg}{Numeric vector contains the ending degree of each bin}
\item{interval_length}{Numeric vector contains the length of each bin.}
\item{binning}{Logical. Indicates whether binning was applied or not.}
\item{g}{Integer. Number of bins}
\item{time_compress_mode}{Integer. The mode of time compression.}
\item{t_compressed}{Integer. The number of time stamps actually used}
\item{compressed_unique_time}{The time stamps that are actually used}
\item{compress_ratio}{Numeric.}
\item{custom_time}{Vector. The time stamps specified by user.}
}

\author{
Thong Pham \email{thongpham@thongpham.net}
}
\seealso{
For creating the needed input for this function (a \code{PAFit_net} object), see \code{\link{as.PAFit_net}}, \code{\link{from_igraph}}, \code{\link{from_networkDynamic}}, and \code{\link{graph_from_file}}. 

For the next step, see \code{\link{Newman}}, \code{\link{Jeong}} or \code{\link{only_A_estimate}} for estimating the attachment function in isolation, \code{\link{only_F_estimate}} for estimating node fitnesses in isolation, and \code{\link{joint_estimate}} for joint estimation of the attachment function and node fitnesses.
}



\examples{
library("PAFit")
net        <- generate_BA(N = 100 , m = 1)
net_stats  <- get_statistics(net)
summary(net_stats)
}

\concept{temporal complex networks}
