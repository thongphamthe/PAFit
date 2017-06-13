graph_from_file <- function(file_name, format = "edgelist", type = "directed") {
  if ("edgelist" == format) {
      graph <- as.matrix(read.table(file_name))  
  } else if ("gml" == format) {
      cat("\nOption type is ignored. Use the type in the gml file.\n")  
      raw_graph <- read_graph(file_name, format = "gml")  
      #convert to the edgelist format of PAFit_net
      graph <- cbind(ends(raw_graph, E(raw_graph)),E(raw_graph)$time)
      if (sum(V(raw_graph)$isolated) > 0) {
          # there is isolated node
          vertex  <- vertex_attr(raw_graph, index = V(raw_graph))
          time    <- vertex$time[vertex$isolated == 1]  
          node_id <- vertex$id[vertex$isolated == 1]
          graph   <- rbind(graph, cbind(node_id,-1,time))
      }
  } else {
      stop("Error: Unrecognized format.")  
  }
  
  result <- list(graph = graph, type = type, PA = NULL, fitness = NULL)
  class(result) <- "PAFit_net"
  return(result)
}