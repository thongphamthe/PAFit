\name{graph_from_file}
\alias{graph_from_file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Read file to a PAFit_net object
}
\description{
  This function reads an input file to a \code{PAFit_net} object. Accepted formats are the edgelist format or the \code{gml} format. 
}
\usage{
 graph_from_file(file_name, format = "edgelist", type = "directed")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_name}{
    A string indicates the file name.
  }
  \item{format}{
    String. Possible values are \code{"edgelist"} and \code{"gml"}.
    
If \code{format} is \code{"edgelist"}, we assume the following edgelist matrix format. Each row is assumed to be of the form (\code{from_node_id} \code{to_node_id} \code{time_stamp}).
\code{from_node_id} is the id of the source node. \code{to_node_id} is the id of the destination node. \code{time_stamp} is the arrival time of the edge. \code{from_node_id} and \code{to_node_id} are assumed to be integers that are at least \eqn{0}. They need not to be contiguous.

To register a new node \eqn{i} at time \eqn{t} without any edge, add a row with format (\code{i -1 t}). This works for both undirected and directed networks.

\code{time_stamp} can be either numeric or string. The value of a time-stamp can be arbitrary, but we assume that a smaller time_stamp (regarded so by the \code{sort} function in \code{R}) represents an earlier arrival time. Examples of time-stamps that satisfy this assumption are the integer \code{0:T}, the string format `yyyy-mm-dd', and the POSIX time.

  If \code{format} is \code{"gml"}, there must be a binary field \code{directed} indicating the type of the network (\code{0}: undirected, \code{1}: directed). The required fields for an edge are: \code{source}, \code{target}, and \code{time}. \code{source} and \code{target} are the ID of the source node and the target node, respectively. \code{time} is the time-stamp of the edge. The required fields for a node are: \code{id}, \code{isolated} (binary) and \code{time}. The binary field \code{isolated} indicates whether this node is an isolated node when it enters the system or not. If \code{isolated} is \code{1}, then \code{time} must contain the node's appearance time. If \code{isolated} is \code{0}, then we can automatically infer the node's appearance time from its edges, so the field \code{time} in this case can be \code{NULL}. The assumptions on node IDs and the format of time-stamps are the same as in the case when \code{format = "edgelist"}. See \code{\link{graph_to_file}} to see detail on the format of the \code{gml} file this package outputs.
}
\item{type}{
   String. Indicates whether the network is \code{"directed"} or \code{"undirected"}. This option is ignored if \code{format} is \code{"gml"}, since the information is assumed to be contained in the \code{gml} file.
  }
}


\value{
  An object of class \code{PAFit_net} containing the network.
}
\author{
  Thong Pham \email{thongpham@thongpham.net}
}


\examples{
  library("PAFit")
  # a network from Bianconi-Barabasi model
  net        <- generate_BB(N = 50 , m = 10 , s = 10)
  
  #graph_to_file(net, file_name = "test.gml", format = "gml")
  #reread    <- graph_from_file(file_name = "test.gml", format = "gml")
}
