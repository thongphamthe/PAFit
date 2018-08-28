plot.PAFit_net <- function(x, 
                           plot = "graph",
                           slice = length(unique(x$graph[,3])) - 1,
                           ...) {
 
  net          <- x
  T            <- length(unique(net$graph[,3]))
  if (slice < 0)
      stop("Error: slice should be non-negative") 
  if (slice + 1 > T)
      cat("\nWarning: slice is bigger than the final time-stamp. Plotted the final snapshot.")
  unique_time  <- sort(unique(net$graph[,3]))
  graph        <- net$graph[order(net$graph[,3]),]
  time         <- graph[,3]
  list_net     <- list(length = T)
  current_node <- NULL 
  if (net$type == "directed") {
    directed <- TRUE
  } else directed <- FALSE
  
  if ("graph" == plot) {
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
  
  
      edge_portion <- graph[time <= unique_time[slice + 1], 1:2, drop = FALSE]
      out_node     <- edge_portion[edge_portion[, 2] != -1, 1] # source node
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
    
      my_graph <- network(ajc_matrix[as.character(as.integer(current_node)), 
                          as.character(as.integer(current_node))], directed = directed, 
                          loops = TRUE, multiple = FALSE)
      plot(my_graph, ...)
  } else if ("degree" == plot) {
      edge_portion <- graph[time <= unique_time[slice + 1], 1:2, drop = FALSE]
      out_node     <- edge_portion[edge_portion[, 2] != -1, 1] # source node
      in_node      <- edge_portion[edge_portion[, 2] != -1, 2] # destination node
      appear_node  <- edge_portion[edge_portion[, 2] == -1, 1] # isolated nodes at this time
      current_node <- sort(unique(c(current_node,union(c(in_node,out_node),appear_node))))
      degree_vec   <- rep(0,length(current_node))
      names(degree_vec) <- as.character(as.integer(current_node))
      if (TRUE == directed) {
          temp       <- table(in_node)  
          degree_vec[as.character(as.integer(labels(temp)[[1]]))]  <- 
                degree_vec[as.character(as.integer(labels(temp)[[1]]))] + temp
             
      } else{
          temp       <- table(c(in_node,out_node))
          degree_vec[as.character(as.integer(labels(temp)[[1]]))]  <- 
            degree_vec[as.character(as.integer(labels(temp)[[1]]))] + temp
      }
      degree_dist <- table(degree_vec)
      degree      <- as.integer(labels(degree_dist)[[1]])
      plot(degree + 1, degree_dist, log = "xy", xlab = "Degree + 1", ylab = "Frequency", pch = 20,
           col = "red", type = "n", axes = FALSE, 
           mgp = c( 2.5 , 1 , 0 ),
           tcl = 0.5,
           ...)
      magaxis(grid = FALSE, frame.plot = TRUE,
              mgp = c( 2.5 , 1 , 0 ),
              tcl = 0.5,...)
      points(degree + 1, degree_dist, pch = 20, col = "red",...)
  } else if ("PA" == plot) {
      # plot the PA function  
      if (is.null(net$PA))
          stop("Error: the object does not contain the PA field.")
    
      deg      <- 0:(length(net$PA) - 1)
      PA       <- net$PA
      non_zero <- which(PA > 10^-20 & deg > 0) 
      
      xlim     <- c(min(deg[non_zero]) , max(deg[non_zero]))  
      ylim     <- c(min(PA[non_zero]), max(PA[non_zero]))  

      col_point   <- "grey25"
      shade_point <- 0.5
      col_pa <- as.vector(col2rgb(col_point)) / 255
    
      plot(deg[non_zero], PA[non_zero], xlab = "Degree k",
           ylab = "Attachment function",
           xlim = xlim , ylim = ylim , 
           axes = FALSE, log = "xy" , 
           col = rgb(col_pa[1],col_pa[2],col_pa[3], shade_point), 
           mgp = c( 2.5 , 1 , 0 ), 
           pch = 20,...)
     magaxis(grid = FALSE, frame.plot = TRUE, usepar=TRUE);
 
    
  } else if ("fit" == plot) {
      # plot the distribution of node fitness
      non_zero <- net$fitness > 10^-20 
    
      f_non <- net$fitness[non_zero]
      d     <- density(f_non)
      ok_d  <- d$x > 0
      max_x <- max(d$x[ok_d])
      min_x <- min(d$x[ok_d])
      
    
      plot(d$x[ok_d] , d$y[ok_d], col  = 2 , log = "x", lwd = 0, 
           main = "", xlab = "Fitness", ylab = "Density",
           mgp = c( 2.5 , 1 , 0 ), 
           axes = FALSE,
           pch = "",...)
      axis(side = 1, at=c(format(round(min_x, 2), nsmall = 2),
                          format(round(1, 2), nsmall = 2),
                          format(round(max_x, 2), nsmall = 2)),
           tcl = 0.5, ...)
      magaxis(side = 2,grid = FALSE, frame.plot = TRUE,
              logpretty = FALSE,...)
      u    <- smooth.spline(d$x, d$y, spar = 0.01)
      ok_u <- u$x > 0 & u$y > 0 
      lines(u$x[ok_u], u$y[ok_u], col = "grey50",lwd = 2.5);
    
  } else {
     stop("Error: Unrecognized plot type.")  
  }
}
