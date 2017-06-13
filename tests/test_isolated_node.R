
rm(list = ls())
library(PAFit)
set.seed(1)
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

stats_new <- get_statistics(net_new)
result_new <- only_A_estimate(net_new, stats_new)

net_approximate <- rbind(net$graph,
                         cbind(isolated_node,isolated_node,random_time))
  
net_stats <- get_statistics(as.PAFit_net(net_approximate)) 
result    <-  only_A_estimate(as.PAFit_net(net_approximate), net_stats)

result     <- result$estimate_result
result_new <- result_new$estimate_result
A_old <- result$theta[result$center_k >= 2]
A_old <- A_old/A_old[1]
A_new <- result_new$theta[result_new$center_k >= 2]
A_new <- A_new / A_new[1]
center_k <- result$center_k[result$center_k >= 2]

plot(center_k, A_old, log = "xy", pch = 20, 
     xlab = "Degree k", col = rgb(1,0,0,0.4), cex = 2,
     ylab = "PA value")

points(center_k,A_new, pch = 16, cex = 2, col = rgb(0,1,0,0.4))
