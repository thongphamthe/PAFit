from_igraph <- function(net) {
  if (!is(net,"igraph")) {
    stop("net must be an object of class igraph.")  
  }
  graph <- matrix(as.integer(as.matrix(as_data_frame(net))), ncol = 3)
  if (is_directed(net))
    type <- "directed"
  else type <- "undirected"
  result <- list(graph = graph, type = type, PA = NULL, fitness = NULL)
  class(result) <- "PAFit_net"
  return(result)
}