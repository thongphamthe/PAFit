to_networkDynamic <- function(net_object) {
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
  }else directed <- FALSE
  
  edge_portion <- graph[, 1:2, drop = FALSE]
  out_node     <- edge_portion[edge_portion[,1] != -1, 1]    # source node
  in_node      <- edge_portion[edge_portion[,2] != -1, 2]    # destination node
  appear_node  <- edge_portion[edge_portion[,2] == -1, 1]    # isolated nodes at this time
  current_node <- unique(c(in_node,out_node,appear_node))
  N            <- length(current_node)
  
  ajc_matrix           <- matrix(0, nrow = N, ncol = N)
  rownames(ajc_matrix) <- as.character(as.integer(current_node))
  colnames(ajc_matrix) <- as.character(as.integer(current_node))
  current_node         <- NULL
 
  for (t in 1:T) {
      edge_portion <- graph[time == unique_time[t], 1:2, drop = FALSE]
      out_node     <- edge_portion[edge_portion[, 1] != -1, 1] # source node
      in_node      <- edge_portion[edge_portion[, 2] != -1, 2] # destination node
      appear_node  <- edge_portion[edge_portion[, 2] == -1, 1] # isolated nodes at this time
      current_node <- sort(unique(c(current_node,union(c(in_node,out_node),appear_node))))
      edge_portion <- as.data.frame(edge_portion)
      repeat_vec   <- ddply(edge_portion,.(edge_portion[,1],edge_portion[,2]), nrow)
      index        <- matrix(c(as.character(repeat_vec[,1]), as.character(repeat_vec[,2])), ncol = 2)
    
      ajc_matrix[index] <- ajc_matrix[index] + repeat_vec[,3]
    
      if (directed == FALSE) {
          index <- matrix(c(as.character(repeat_vec[,2]), as.character(repeat_vec[,1])), ncol = 2)  
          ajc_matrix[index] <- ajc_matrix[index] + repeat_vec[,3]
      }
          
      list_net[[t]]  <- as.network(ajc_matrix[as.character(as.integer(current_node)), 
                                              as.character(as.integer(current_node))], 
                                   directed = directed, loops = TRUE, multiple = TRUE)
  }
  # the whole to_networkDynamic function is slow, but the bottleneck is in the final step: converting list of nets to networkDynamic
  # so there is little things I can do now
  result <- networkDynamic(network.list = list_net, vertex.pid = "vertex.names", verbose = FALSE)
  return(result)
}

