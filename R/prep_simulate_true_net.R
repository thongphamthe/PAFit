
.prep_simulate_true_net <- function(net, net_type = "directed"){
  net               <- net[order(net[,3], decreasing = FALSE),]
  time_stamp        <- as.vector(net[,3])
  in_node           <- as.vector(net[,2])
  out_node          <- as.vector(net[,1])
  out_node          <- out_node
  node_id           <- as.integer(sort(union(in_node[in_node !=  -1],out_node[out_node != - 1])))
  
  unique_time       <- sort(unique(time_stamp))
  T                 <- length(unique_time)
  N                 <- length(node_id)
  max_node_id       <- max(node_id);

  appear            <- rep(0,N)
  names(appear)     <- as.integer(node_id)
  
  m_t               <- rep(0, T - 1);
  isolated          <- list(length = T - 1);
  start_point       <- list(length = T - 1); # this will be empty in the case of undirected network
  new_but_not_count <- list(length = T - 1); 
  
  G_0                          <- as.matrix(net[time_stamp == unique_time[1],,drop = FALSE])
  node_0                       <- as.vector(G_0[,1:2])
  node_0                       <- node_0[node_0 != -1]
  appear[as.character(as.integer(node_0))] <- 1
  
  for (t in 2:T) {
      G_t     <- net[time_stamp == unique_time[t], ,drop = FALSE]
      #node_t  <- as.vector(G_t[,1:2])
      out_node_t <- G_t[,1]
      in_node_t  <- G_t[,2]
      isolated[[t - 1]] <- out_node_t[which(in_node_t == -1)]
      
      out_node_t    <- out_node_t[which(in_node_t != -1)]
      
      G_t           <- G_t[which(in_node_t != -1),,drop = FALSE]
      in_node_t     <- in_node_t[which(in_node_t != -1)]
      
      if (net_type == "directed")
          new_but_not_count[[t - 1]] <- G_t[appear[as.character(as.integer(in_node_t))] == 0,,drop = FALSE]
      else {
          new_but_not_count[[t - 1]] <- G_t[appear[as.character(as.integer(in_node_t))] == 0 | 
                                            appear[as.character(as.integer(out_node_t))] == 0, , drop = FALSE]   
      }
      if (net_type == "directed") {
           m_t[t - 1] <- sum(appear[as.character(as.integer(in_node_t))] == 1)
      } else   m_t[t - 1] <- sum(appear[as.character(as.integer(in_node_t))] == 1 & appear[as.character(as.integer(out_node_t))] == 1)
      
      if (net_type == "directed")
          start_point[[t - 1]]       <- out_node_t[appear[as.character(as.integer(in_node_t))] == 1]
      
      appear[as.character(as.integer(in_node_t))]  <- 1
      appear[as.character(as.integer(out_node_t))] <- 1
  }
  result <- list(G_0 = G_0, m_t = m_t, start_point = start_point, new_but_not_count = new_but_not_count,
                 isolated = isolated, node_id = node_id, net_type = net_type,
                 unique_time = unique_time)
  return(result)
}