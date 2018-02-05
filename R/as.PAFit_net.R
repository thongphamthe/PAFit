
as.PAFit_net <- function(graph, type = "directed", PA = NULL, fitness = NULL) {
 if (dim(graph)[2] != 3) 
     stop("Error: graph should have three columns.")
 if (type != "directed" && type != "undirected")
      stop("Error: type should be either directed or undirected.")
 if (!is.null(PA)){
    if (sum(is.na(PA)) > 0)
        stop("Error: There is NAs in the PA function.")        
    if (sum(PA <= 0)  > 0)
        stop("Error: PA function should be positive.")    
 }
 if (!is.null(fitness)){
    if (sum(is.na(fitness)) > 0)
      stop("Error: There is NAs in node fitnesses.")        
    if (sum(fitness <= 0)  > 0)
      stop("Error: Node fitnesses should be positive.")    
 }
 ### sort the edges based on their arrival time
 graph  <- graph[order(graph[,3]),]
 result <- list(graph = graph, type = type, PA = PA, fitness = fitness)
 class(result) <- "PAFit_net"
 return(result)
}
