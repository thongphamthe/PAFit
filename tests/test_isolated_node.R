
rm(list = ls())
library(PAFit)

# for CRAN. In developing, set ii from 1 to 1000
for (ii in 1) {
    set.seed(2)  
    print(ii)
    prob_m    <- "FALSE"
    inc       <- "FALSE"
    log       <-  c("FALSE")
    net       <-  generate_BA(N = 10, m = 5)


    max_id        <- max(net$graph)
    random_time   <- sample(1:max(net$graph[,3]),size = 50, replace = TRUE)
    isolated_node <- (max_id + 1):(max_id+1 + 50 - 1) 

    net_new <- net
    net_new$graph <- rbind(net$graph,
                       cbind(isolated_node,rep(-1,50),random_time))

    stats_new <- get_statistics(net_new, g = 3)
    result_new <- only_A_estimate(net_new, stats_new, stop_cond = 10^-2)
    result_temp_temp <- joint_estimate(net_new, stats_new, stop_cond = 10^-2)
    
    stats_new <- get_statistics(net_new)
    result_new <- only_A_estimate(net_new, stats_new, stop_cond = 10^-2)
  

    net_approximate <- rbind(net$graph,
                         cbind(isolated_node,isolated_node,random_time))
  
    net_stats <- get_statistics(as.PAFit_net(net_approximate)) 
    result    <-  only_A_estimate(as.PAFit_net(net_approximate), net_stats, stop_cond = 10^-2)
    result_temp <- joint_estimate(as.PAFit_net(net_approximate), net_stats, stop_cond = 10^-2)

    result     <- result$estimate_result
    result_new <- result_new$estimate_result
    A_old <- result$theta[result$center_k >= 2]
   A_old <- A_old/A_old[1]
   A_new <- result_new$theta[result_new$center_k >= 2]
    A_new <- A_new / A_new[1]
    center_k <- result$center_k[result$center_k >= 2]
   #plot(center_k, A_old, log = "xy", pch = 20, 
  #   xlab = "Degree k", col = rgb(1,0,0,0.4), cex = 2,
  #   ylab = "PA value")
   # points(center_k,A_new, pch = 16, cex = 2, col = rgb(0,1,0,0.4))
}

