print.PAFit_net <- function(x,...) {
 cat(paste0("An object of class PAFit_net containing a dynamic network"))
  cat(paste0("\nType: ",x$type)) 
 cat(paste0("\nNumber of unique time points: ", length(unique(x$graph[,3]))))
 nodes <- unique(c(as.vector(x$graph[,1]), as.vector(x$graph[x$graph[,2] != -1,2])))
 num_edges <- sum(x$graph[,2] != -1)
 cat(paste0("\nNumber of nodes: ", length(nodes)))
 cat(paste0("\nNumber of edges: ", num_edges))
 cat("\n")
}