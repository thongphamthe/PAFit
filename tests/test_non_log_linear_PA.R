#testing CV-related functions

## code here that contains the tests


rm(list = ls())
library(PAFit)
for (ii in 1) {
  set.seed(ii)
  print(ii)
  prob_m <- "FALSE"
  inc    <- "FALSE"
  log    <-  c("FALSE")
  beta <- 1
  alpha <- 2
  net  <- generate_net(N = 20, m = 5,prob_m = prob_m, 
                       increase = inc, log = log, multiple_node = 2, num_seed = 10,
                       mode = 3, s = 10,alpha = alpha, beta = beta)
  
  net_stats <- get_statistics(net,deg_threshold = 0, 
                              binning = TRUE, g = 50) 
  
  result <- joint_estimate(net, net_stats)
  k      <- result$estimate_result$k
  true_A <- alpha * log(pmax(1,k))^beta + 1
  estimate_A <- result$estimate_result$A
  estimate_A <- estimate_A/estimate_A[1]
  
  max_A  <- max(true_A, result$estimate_result$A)
  
  #plot(result,net_stats)

  plot(k + 1,estimate_A, ylim = c(1,max_A),
       log = "xy")
  lines(k + 1,true_A)
}




