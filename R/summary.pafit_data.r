summary.PAFit_data <- function(object,...){
  print(object)
  cat("Type of network: ",object$net_type[1],"\n");
  cat("Number of nodes in the final network: ",length(object$final_deg),"\n")
  cat("Number of edges in the final network: ",sum(object$final_deg),"\n")
  cat("Number of new nodes: ",length(object$final_deg) - object$initial_nodes,"\n")
  cat("Number of new edges: ",sum(object$sum_m_k), "\n")
  cat("Numbef of time-steps:",  object$T,"\n");
  cat("Maximum degree: ", object$deg_max,"\n");
  cat("Number of bins:", object$g,"\n");
}