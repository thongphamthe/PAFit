to_igraph <- function(net_object) {
  net <- net_object
  if (!is(net,"PAFit_net")) {
    stop("net must be an object of class PAFit_net.")  
  }
  T            <- length(unique(net$graph[,3]))
  unique_time  <- sort(unique(net$graph[,3]))
  graph        <- net$graph[order(net$graph[,3]),]
  time         <- graph[,3]
  list_net     <- list(length = T)
  current_node <- NULL 
  if (net$type == "directed") {
      directed <- TRUE
  } else directed <- FALSE
  colnames(net$graph) <- c("id_1","id_2","time")
  result <- graph.data.frame(net$graph, directed = directed)
  return(result)
}

