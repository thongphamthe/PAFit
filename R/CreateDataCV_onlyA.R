
.CreateDataCV_onlyA <- function(net                   , p          = 0.75 , g           = 50, 
                               net_type = "directed" , deg_thresh = 0    ) {
  net               <- net[order(net[,3] , decreasing = FALSE),]
  time_stamp        <- as.vector(net[,3])
  in_node           <- as.vector(net[,2])
  
  out_node          <- as.vector(net[,1])
  node_id           <- sort(union(in_node,out_node))
  
  names(node_id)    <- node_id
  unique_time       <- sort(unique(time_stamp))
  
  T                 <- length(unique_time)
  
  N                 <- length(node_id) 
    
  first_time        <- time_stamp[1]
  edge_cumsum       <- cumsum(as.vector(table(time_stamp[time_stamp != first_time]))) 
  edge_ratio        <- edge_cumsum/edge_cumsum[length(edge_cumsum)]
  use_time          <- unique_time[which(edge_ratio >= p)[1]]
  
  data_new          <- net[time_stamp <= use_time, ]
  stats             <- get_statistics(data_new , net_type = net_type, 
                                     binning = TRUE , g = g , deg_threshold = deg_thresh)
  
  final_stat        <- get_statistics(net, binning = TRUE, g = g, net_type = net_type, deg_threshold = 0)
  
  n_tk_each         <- final_stat$n_tk[unique_time[unique_time > use_time & unique_time < T],]
  
  m_each            <- final_stat$m_t[unique_time[unique_time > use_time & unique_time < T]]
  
  m_tk_each         <- final_stat$m_tk[unique_time[unique_time > use_time & unique_time < T],]
  
  prob_em_each      <- m_tk_each / m_each
 
 
  result            <- list(stats        = stats              , n_tk_each = n_tk_each ,
                            m_each       = m_each             , p         = p         ,
                            center_k     = final_stat$center_k,  
                            prob_em_each = prob_em_each       , use_time  = use_time)
  class(result)     <- "CV_Data"

  return(result)
}


