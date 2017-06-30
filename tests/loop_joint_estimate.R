

if (FALSE) {
  #setwd("tests")
  rm(list = ls())
  library(PAFit)
  set.seed(4)
  prob_m <- "FALSE"
  inc    <- "FALSE"
  log    <-  c("FALSE")
  
  M <- 10
  
  alpha_vec     <- rep(0,M)
  s_vec         <- rep(0,M)
  alpha_optimal <- rep(0,M)
  r_optimal     <- rep(0,M)
  for (i in 1:M) {
      #set.seed(1)
      net  <- generate_BB(N = 1000, m = 15, num_seed = 200, 
                          multiple_node = 200,
                          s = 1) 
                          #shape = 1, rate = 1)
  
      net_stats <- get_statistics(net, deg_threshold = 1) 
  
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
