\name{graph_to_file}
\alias{graph_to_file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Write the graph in a PAFit_net object to file
}
\description{
  This function writes a graph in a \code{PAFit_net} object to an output file. Accepted file formats are the edgelist format or the \code{gml} format. 
  
}
\usage{
graph_to_file(net_object, file_name, format = "edgelist")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{net_object}{
An object of class \code{PAFit_net}.
}
\item{file_name}{
A string indicates the file name.
}
\item{format}{
String. Possible values are \code{"edgelist"} and \code{"gml"}.

If \code{format = "edgelist"}, we just output the edgelist matrix contained in the \code{PAFit_net} object as it is.

If \code{format = "gml"}, here is the specification of the \code{gml} file. There is a binary field \code{directed} indicating the type of the network (\code{0}: undirected, \code{1}: directed). There are three atrributes for an edge: \code{source}, \code{target}, and \code{time}. There are three atrributes for a node: \code{id}, \code{isolated} (binary) and \code{time}. The atrribute \code{time} is \code{NULL} if the attribute \code{isolated} is \code{0} (since this is not an isolated node, we do not need to record its first apperance time). On the other hand, \code{time} is the node's appearance time if attribute \code{isolated} is \code{1}.
}
}

\value{
The function writes directly to the output file.
}

\author{
Thong Pham \email{thongpham@thongpham.net}
}


\examples{
library("PAFit")
# a network from Bianconi-Barabasi model
net        <- generate_BB(N = 50 , m = 10 , s = 10)
#graph_to_file(net, file_name = "test.gml", format = "gml")
}
