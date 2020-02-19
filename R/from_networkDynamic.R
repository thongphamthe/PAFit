
from_networkDynamic <- function(net) {
  if (!is(net,"networkDynamic")) {
    stop("net must be an object of networkDynamic.")  
  }
  result        <- as.data.frame(net)
  graph         <- cbind(result$tail,result$head,result$onset)
  if (is.directed(net) == TRUE) {
      type <- "directed"
  } else type <- "undirected"
  return(as.PAFit_net(graph = graph, type = type))
}

# vertex.attr = list("name" = as.character(as.integer(current_node))), 
# vertex.attrnames = list("name")