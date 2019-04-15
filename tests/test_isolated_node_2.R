
rm(list = ls())
library(PAFit)

# for CRAN. In developing, set ii from 1 to 1000
for (ii in 1) {
  set.seed(2) 
  print(ii)
  prob_m    <- "FALSE"
  inc       <- "FALSE"
  log       <-  c("FALSE")
  net       <-  generate_BA(N = 100, m = 5)


  max_id        <- max(net$graph)
  random_time   <- sample(1:max(net$graph[,3]),size = 500, replace = TRUE)
  isolated_node <- (max_id + 1):(max_id+1 + 500 - 1) 

  net_new <- net
  net_new$graph <- rbind(net$graph,
                       cbind(isolated_node,rep(-1,500),random_time))

  new_edge <- sample(1:max_id, size = 500, replace = TRUE)
  T <- max(net$graph[,3])

  net_new$graph <- rbind(net_new$graph, cbind(isolated_node,new_edge, T + 1))

  stats_new <- get_statistics(net_new)
  result_new <- only_A_estimate(net_new,stats_new, stop_cond = 10^-3)
  result <- joint_estimate(net_new, stats_new, stop_cond = 10^-3)

  #plot(result_new,stats_new)
}

