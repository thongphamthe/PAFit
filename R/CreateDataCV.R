
.CreateDataCV<- function(net_object               , p          = 0.75 , g           = 50, 
                        net_type = "directed" , deg_thresh = 0    , exclude_end = FALSE) {
  #net               <- as.matrix(net)
  oopts <- options(scipen = 999)
  on.exit(options(oopts))
  net               <- net_object$graph
  net_type          <- net_object$type
  net               <- net[order(net[,3] , decreasing = FALSE),]
  time_stamp        <- as.vector(net[,3])
  time_stamp        <- time_stamp[net[,2] != -1]
  in_node           <- as.vector(net[,2])
  in_node           <- in_node[net[,2] != -1]
  out_node          <- as.vector(net[,1])
  out_node          <- out_node[net[,2] != -1]
  
  ### use as.numeric, instead of as.integer, since the node id might be outside of integer range
  node_id           <- as.numeric(sort(union(in_node,out_node)))
  names(node_id)    <- as.numeric(node_id)
  unique_time       <- sort(unique(time_stamp))
  T                 <- length(unique_time)
  N                 <- length(node_id) 
  appear            <- rep(0 , N)
  names(appear)     <- as.numeric(node_id)
  first_time        <- unique_time[1]
  edge_cumsum       <- cumsum(as.vector(table(time_stamp[time_stamp != first_time]))) 
  edge_ratio        <- edge_cumsum/edge_cumsum[length(edge_cumsum)]
  
  ok_time           <- which(edge_ratio >= p)
  
  if (length(ok_time) == 1) {  ## must check, since base on p, all the data might be regarded as learning data
    use_time      <- unique_time[length(unique_time) - 1]
  } else
    use_time      <- unique_time[which(edge_ratio >= p)[1]]
  
  data_new          <- net[time_stamp <= use_time, ]
  net_new           <- as.PAFit_net(graph = data_new, type = net_type)
  stats             <- get_statistics(net_new,
                                      binning = TRUE , g = g , deg_threshold = deg_thresh)
  appear[as.character(as.numeric(stats$f_position))] <- 1
  deg                                    <- stats$final_deg[as.character(as.numeric(stats$f_position))]
  if (FALSE == exclude_end) {
      if (net_type == "directed") {
          prob_em_each              <- matrix(0 , nrow = sum(unique_time > use_time) , ncol = length(stats$f_position))
          colnames(prob_em_each)    <- as.numeric(stats$f_position)
          m_each                    <- rep(0 , sum(unique_time > use_time))
          deg_each                  <- matrix(0 , nrow = sum(unique_time > use_time) , ncol = length(stats$f_position))
          colnames(deg_each)        <- as.numeric(stats$f_position)
          deg_each[1,]              <- deg
          time_each                 <- unique_time[unique_time > use_time]
          for (i in 1:length(time_each)){
              new_links                            <- in_node[time_stamp == time_each[i]]
              new_links                            <- new_links[appear[as.character(as.numeric(new_links))] == 1]
              m_each[i]                            <- length(new_links)
              aaa                                  <- table(new_links)
              prob_em_each[i , as.character(as.numeric(labels(aaa)[[1]]))]   <- aaa
              prob_em_each[i , as.character(as.numeric(labels(aaa)[[1]]))]   <- prob_em_each[i, as.character(as.numeric(labels(aaa)[[1]]))]/ m_each[i] 
          if (i < length(time_each)) {
              deg_each[i + 1,]                     <- deg_each[i,];    
              deg_each[i + 1,as.character(as.numeric(labels(aaa)[[1]]))]     <- deg_each[i + 1 , as.character(as.numeric(labels(aaa)[[1]]))] + aaa
          }
      }
  } else { #undirected network
        prob_em_each              <- matrix(0,nrow = sum(unique_time > use_time),ncol = length(stats$f_position))
        colnames(prob_em_each)    <- as.numeric(stats$f_position)
        m_each                    <- rep(0,sum(unique_time > use_time))
        deg_each                  <- matrix(0,nrow = sum(unique_time > use_time),ncol = length(stats$f_position))
        colnames(deg_each)        <- as.numeric(stats$f_position)
        deg_each[1,]              <- deg
        time_each                 <- unique_time[unique_time > use_time]
        for (i in 1:length(time_each)){
            new_in_links      <- in_node[time_stamp == time_each[i]]
            new_in_links      <- new_in_links[appear[as.character(as.numeric(new_in_links))] == 1]
            new_out_links     <- out_node[time_stamp == time_each[i]]
            new_out_links     <- new_out_links[appear[as.character(as.numeric(new_out_links))] == 1]
            m_each[i]         <- length(c(new_in_links,new_out_links))
            aaa               <- table(c(new_in_links,new_out_links))
            prob_em_each[i, as.character(as.numeric(labels(aaa)[[1]]))] <- aaa
            prob_em_each[i, as.character(as.numeric(labels(aaa)[[1]]))] <- prob_em_each[i, as.character(as.numeric(labels(aaa)[[1]]))]/ m_each[i] 
            if (i < length(time_each)) {
                deg_each[i+1,]                 <- deg_each[i,];    
                deg_each[i+1, as.character(as.numeric(labels(aaa)[[1]]))] <- deg_each[i+1, as.character(as.numeric(labels(aaa)[[1]]))] + aaa
        }
      }
    
  }
  result            <- list(stats        = stats       , deg_each = deg_each,
                            m_each       = m_each      , p        = p       ,
                            prob_em_each = prob_em_each, use_time = use_time)
  class(result)     <- "CV_Data"
  return(result)
  } else {
  #exclude_end = TRUE
      deg_max <- stats$deg.max
      if (net_type == "directed") {
          prob_em_each              <- matrix(-1 , nrow = sum(unique_time > use_time) , ncol = length(stats$f_position))
          colnames(prob_em_each)    <- as.numeric(stats$f_position)
          m_each                    <- rep(-1 , sum(unique_time > use_time))
          deg_each                  <- matrix(-1 , nrow = sum(unique_time > use_time) , ncol = length(stats$f_position))
          colnames(deg_each)        <- as.numeric(stats$f_position)
      
 
          deg_each[1,]              <- deg
          time_each                 <- unique_time[unique_time > use_time]
          for (i in 1:length(time_each)){
              new_links       <- in_node[time_stamp == time_each[i]]
              new_links       <- new_links[appear[as.character(as.numeric(new_links))] == 1]
              ### IMPORTANT: remove nodes that go outside the range of the degree distribution of testing net
              new_links_final <- new_links[deg_each[i,as.character(as.numeric(new_links))] < deg_max]
              m_each[i]       <- length(new_links_final)
              aaa             <- table(new_links)
              bbb             <- table(new_links_final)
              prob_em_each[i, as.character(as.numeric(labels(bbb)[[1]]))] <- bbb
              prob_em_each[i, as.character(as.numeric(labels(bbb)[[1]]))] <- prob_em_each[i, as.character(as.numeric(labels(bbb)[[1]])) ]/ m_each[i] 
              
              if (i < length(time_each)) {
                  deg_each[i + 1 , ]                     <- deg_each[i , ];    
                  deg_each[i + 1 , as.character(as.numeric(labels(aaa)[[1]]))]     <- deg_each[i + 1 , as.character(as.numeric(labels(aaa)[[1]]))] + aaa
              }
      }
      result            <- list(stats        = stats       , deg_each = deg_each,
                                m_each       = m_each      , p        = p       ,
                                prob_em_each = prob_em_each, use_time = use_time)
      class(result)     <- "CV_Data"
      return(result)
    }
  }
}


