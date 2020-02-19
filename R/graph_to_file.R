graph_to_file <- function(net_object, file_name, format = "edgelist") {
  net <- net_object
  #type = "edgelist" or "gml"
  if (!is(net_object,"PAFit_net")) {
      stop("Error: net must be an object of class PAFit_net.")  
  }
  
  if ("edgelist" == format) {
    write.table(net$graph, file = file_name, col.names = FALSE, row.names = FALSE)  
  } else if ("gml" == format) {
    fileConn <- file_name
    cat("graph [\n", file = fileConn, append = FALSE)
    
    if ("directed" == net$type )
        cat("directed 1\n", file = fileConn, append = TRUE)
    else cat("directed 0\n", file = fileConn, append = TRUE)
    
    # write list of nodes
    
    # first get the list of nodes and their appearances time
    
    all_nodes          <- sort(unique(c(as.vector(net$graph[net$graph[,2] != -1,1:2]),
                                      as.vector(net$graph[net$graph[,2] == -1, 1]))))
    isolated_node      <- net$graph[net$graph[,2] == -1, 1]
    time_isolated_node <- net$graph[net$graph[,2] == -1, 3]
  
    non_isolated_node  <- setdiff(all_nodes,isolated_node)
    if (length(non_isolated_node) > 0)
    for (v in 1:length(non_isolated_node)) {
        cat("node [\n", file = fileConn, append = TRUE)  
        cat(paste0("id ",non_isolated_node[v],"\n"), file = fileConn, append = TRUE)
        cat("time NULL\n", file = fileConn, append = TRUE)
        cat("isolated 0\n", file = fileConn, append = TRUE)
        cat("]\n", file = fileConn, append = TRUE)
    }
    if (length(isolated_node) > 0)
    for (v in 1:length(isolated_node)) {
      cat("node [\n", file = fileConn , append = TRUE)  
      cat(paste0("id ",isolated_node[v],"\n"), file = fileConn, append = TRUE)
      cat(paste0("time ",time_isolated_node[v],"\n"), file = fileConn, append = TRUE)
      cat("isolated 1\n", file = fileConn, append = TRUE)
      cat("]\n", file = fileConn, append = TRUE)
    }
    
    # write list of edges
    if (dim(net$graph)[1] > 0)
    for (e in 1:dim(net$graph)[1]) 
        if (net$graph[e,2] != -1) { 
        # this is not an isolated node
            cat("edge [\n", file = fileConn, append = TRUE)
            cat(paste0("source ", net$graph[e,1],"\n"), file = fileConn, append = TRUE)
            cat(paste0("target ", net$graph[e,2],"\n"), file = fileConn, append = TRUE)
            cat(paste0("time "  , net$graph[e,3],"\n"), file = fileConn, append = TRUE)
            cat("]\n", file = fileConn, append = TRUE)
        }
    
    
    cat("]", file = fileConn, append = TRUE)
  
  } else {
    stop("Error: unrecognized format")  
  
  } 
}