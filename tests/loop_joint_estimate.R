

if (FALSE) {
  #setwd("tests")
  rm(list = ls())
  library(PAFit)
  set.seed(1)
  prob_m <- "FALSE"
  inc    <- "FALSE"
  log    <-  c("FALSE")
  
  M <- 5
  
  alpha_vec     <- rep(0,M)
  s_vec         <- rep(0,M)
  alpha_optimal <- rep(0,M)
  r_optimal     <- rep(0,M)
  for (i in 1:M) {
      #set.seed(1)
      net  <- generate_net(N = 1000, m = 50, num_seed = 500, 
                          multiple_node = 50, alpha = 1,
                          s = 5) 
                          #shape = 1, rate = 1)
  
      net_stats    <- get_statistics(net) 
  
      print(result <- joint_estimate(net, net_stats))
      
      #result
                      
      alpha_vec[i]     <- result$estimate_result$alpha
      s_vec[i]         <- result$estimate_result$shape
      alpha_optimal[i] <- result$cv_result$alpha_optimal
      r_optimal[i]     <- result$cv_result$r_optimal
  }
  print(alpha_vec)  
  print(mean(alpha_vec))
  print(s_vec)
  print(mean(s_vec))
}
