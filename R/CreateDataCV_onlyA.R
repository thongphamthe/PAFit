
.CreateDataCV_onlyA <- function(net_object , p = 0.75 , g = 50 , deg_thresh = 0    ) {
  #net               <- as.matrix(net)
  oopts <- options(scipen = 999)
  on.exit(options(oopts))
  net               <- net_object$graph
  net_type          <- net_object$type
  net               <- net[order(net[,3] , decreasing = FALSE),]
  time_stamp        <- as.vector(net[,3])
  in_node           <- as.vector(net[,2])
  
  out_node          <- as.vector(net[,1])
  node_id           <- as.numeric(sort(union(in_node,out_node)))
  
  names(node_id)    <- as.numeric(node_id)
  unique_time       <- sort(unique(time_stamp))
  
  T                 <- length(unique_time)
  
  N                 <- length(node_id) 
    
  first_time        <- unique_time[1]
  edge_cumsum       <- cumsum(as.vector(table(time_stamp[time_stamp != first_time]))) 
  edge_ratio        <- edge_cumsum/edge_cumsum[length(edge_cumsum)]
  ok_time           <- which(edge_ratio >= p)
  if (length(ok_time) == 1) {
      use_time      <- unique_time[length(unique_time) - 1]
  } else
      use_time      <- unique_time[which(edge_ratio >= p)[1] + 1]
  
  data_new          <- net[time_stamp <= use_time, ]
  net_new           <- as.PAFit_net(graph = data_new, type = net_type)
  stats             <- get_statistics(net_new,
                                      binning = TRUE , g = g , deg_threshold = deg_thresh)
  
  final_stat        <- get_statistics(net_object, binning = TRUE, g = g, deg_threshold = 0)
  
  n_tk_each         <- final_stat$n_tk[unique_time[-1] > use_time,,drop = FALSE]
  
  m_each            <- final_stat$m_t[unique_time[-1] > use_time]
  
  m_tk_each         <- final_stat$m_tk[unique_time[-1] > use_time,, drop = FALSE]
  
  prob_em_each      <- m_tk_each / m_each
 
 
  result            <- list(stats        = stats              , n_tk_each = n_tk_each ,
                            m_each       = m_each             , p         = p         ,
                            center_k     = final_stat$center_k,  
                            prob_em_each = prob_em_each       , use_time  = use_time)
  class(result)     <- "CV_Data"

  return(result)
}


