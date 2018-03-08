library(PAFit)

net  <- generate_net(N = 50, m = 10,prob_m = TRUE, increase = TRUE, log = TRUE,
                    mode = 1, s = 100)

for (bin in c("FALSE","TRUE")) {
  deg_thresh <- 0
  net_stats <- get_statistics(net,deg_threshold = deg_thresh, binning = bin, g = 10, compress_mode = 2)
  net_stats <- get_statistics(net,deg_threshold = deg_thresh, binning = bin, g = 10, compress_mode = 3, 
                              custom_time = 0:20) 
  
  net_stats <- get_statistics(net,deg_threshold = deg_thresh, binning = bin, g = 10, compress_mode = 3, custom_time = 0:20)
  #check stats
  
  if (sum(net_stats$m_t) != sum(net_stats$sum_m_k))
    stop("wrong at m_t and sum_m_k")
  if (sum(abs(colSums(net_stats$m_tk) - net_stats$sum_m_k)) != 0)
    stop("wrong at m_tk and sum_m_k")
  temp <- sapply(1:(net_stats$t_compressed - 1),function(x) sum(net_stats$node_degree[x,] != -1))
  if (sum(abs(rowSums(net_stats$n_tk) - (rowSums(net_stats$offset_tk) + temp))))
    stop("wrong at node_degree, n_tk, offset_tk") 
  if (sum(net_stats$z_j) > sum(net_stats$m_t))
    stop("wrong at z_j")   
  if (sum(net_stats$z_j) + sum(net_stats$offset_m_tk) - sum(net_stats$sum_m_k))
    stop("Wrong at offset_m_tk") 
}




