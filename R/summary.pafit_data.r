summary.PAFit_data <- function(object,...){
  print(object)
  cat("Type of network:",object$net_type[1],"\n");
  cat("Number of nodes in the final network:",length(object$final_deg),"\n")
  if (object$net_type[1] == "directed") {
      cat("Number of edges in the final network:",sum(object$final_deg),"\n")
  } else {
      cat("Number of edges in the final network:",sum(object$final_deg)/2,"\n")
  }
  cat("Number of new nodes:",length(object$final_deg) - object$initial_nodes,"\n")
  if (object$net_type[1] == "directed") {
    cat("Number of new edges:", sum(object$sum_m_k) + sum(object$offset_m_tk),"\n")
  } else {
      cat("Number of new edges:", sum(object$sum_m_k) + sum(object$offset_m_tk), "\n")
  }
  cat("Number of time-steps:",  object$T,"\n");
  if ("directed" == object$net_type[1]) {
      cat("Maximum in-degree:", object$deg_max,"\n");
  } else {
      cat("Maximum degree:", object$deg_max,"\n");  
  }
  cat("Number of bins:", object$g,"\n");
}