library(PAFit)

net  <- GenerateNet(N = 50, m = 10,prob_m = TRUE, increase = TRUE, log = TRUE,
                    mode = 1, shape = 100, rate = 100)

for (bin in c("FALSE","TRUE")) {
    deg_thresh <- 0
    net_stats <- GetStatistics(net$graph,deg_threshold = deg_thresh, Binning = bin, G = 10, net_type = "undirected") 
            #check stats
    if (sum(net_stats$m_t) != sum(net_stats$Sum_m_k))
        stop("wrong at m_t and sum_m_k")
    if (sum(abs(colSums(net_stats$m_tk) - net_stats$Sum_m_k)) != 0)
        stop("wrong at m_tk and sum_m_k")
    temp <- sapply(1:(net_stats$T-1),function(x) sum(net_stats$node_degree[x,] != -1))
    if (sum(abs(rowSums(net_stats$n_tk) - (rowSums(net_stats$offset_tk) + temp))))
        stop("wrong at node_degree, n_tk, offset_tk") 
    if (sum(net_stats$z_j) > sum(net_stats$m_t))
        stop("wrong at z_j")   
    if (sum(net_stats$z_j) + sum(net_stats$offset_m_tk) - sum(net_stats$Sum_m_k))
        stop("Wrong at offset_m_tk") 
}

net  <- GenerateNet(N = 50, m = 10,prob_m = TRUE, increase = TRUE, log = TRUE,
                    mode = 1, shape = 100, rate = 100)

temp <- net$graph[,2]
net$graph[,2] <- net$graph[,1]
net$graph[,1] <- temp

for (bin in c("FALSE","TRUE")) {
  deg_thresh <- 0
  net_stats <- GetStatistics(net$graph,deg_threshold = deg_thresh, Binning = bin, G = 10, net_type = "undirected") 
  #check stats
  if (sum(net_stats$m_t) != sum(net_stats$Sum_m_k))
    stop("wrong at m_t and sum_m_k")
  if (sum(abs(colSums(net_stats$m_tk) - net_stats$Sum_m_k)) != 0)
    stop("wrong at m_tk and sum_m_k")
  temp <- sapply(1:(net_stats$T-1),function(x) sum(net_stats$node_degree[x,] != -1))
  if (sum(abs(rowSums(net_stats$n_tk) - (rowSums(net_stats$offset_tk) + temp))))
    stop("wrong at node_degree, n_tk, offset_tk") 
  if (sum(net_stats$z_j) > sum(net_stats$m_t))
    stop("wrong at z_j")   
  if (sum(net_stats$z_j) + sum(net_stats$offset_m_tk) - sum(net_stats$Sum_m_k))
    stop("Wrong at offset_m_tk") 
}

net$graph[,1] <- rev(net$graph[,1])

for (bin in c("FALSE","TRUE")) {
  deg_thresh <- 0
  net_stats <- GetStatistics(net$graph,deg_threshold = deg_thresh, Binning = bin, G = 10, net_type = "undirected") 
  #check stats
  if (sum(net_stats$m_t) != sum(net_stats$Sum_m_k))
    stop("wrong at m_t and sum_m_k")
  if (sum(abs(colSums(net_stats$m_tk) - net_stats$Sum_m_k)) != 0)
    stop("wrong at m_tk and sum_m_k")
  temp <- sapply(1:(net_stats$T-1),function(x) sum(net_stats$node_degree[x,] != -1))
  if (sum(abs(rowSums(net_stats$n_tk) - (rowSums(net_stats$offset_tk) + temp))))
    stop("wrong at node_degree, n_tk, offset_tk") 
  if (sum(net_stats$z_j) > sum(net_stats$m_t))
    stop("wrong at z_j")   
  if (sum(net_stats$z_j) + sum(net_stats$offset_m_tk) - sum(net_stats$Sum_m_k))
    stop("Wrong at offset_m_tk") 
}
net$graph[,1] <- rev(net$graph[,1])
net$graph[,2] <- rev(net$graph[,2])

for (bin in c("FALSE","TRUE")) {
  deg_thresh <- 0
  net_stats <- GetStatistics(net$graph,deg_threshold = deg_thresh, Binning = bin, G = 10, net_type = "undirected") 
  #check stats
  if (sum(net_stats$m_t) != sum(net_stats$Sum_m_k))
    stop("wrong at m_t and sum_m_k")
  if (sum(abs(colSums(net_stats$m_tk) - net_stats$Sum_m_k)) != 0)
    stop("wrong at m_tk and sum_m_k")
  temp <- sapply(1:(net_stats$T-1),function(x) sum(net_stats$node_degree[x,] != -1))
  if (sum(abs(rowSums(net_stats$n_tk) - (rowSums(net_stats$offset_tk) + temp))))
    stop("wrong at node_degree, n_tk, offset_tk") 
  if (sum(net_stats$z_j) > sum(net_stats$m_t))
    stop("wrong at z_j")   
  if (sum(net_stats$z_j) + sum(net_stats$offset_m_tk) - sum(net_stats$Sum_m_k))
    stop("Wrong at offset_m_tk") 
}



