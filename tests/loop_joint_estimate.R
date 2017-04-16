if (FALSE) {
  #setwd("tests")
  rm(list = ls())
  library(PAFit)
  set.seed(1)
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
      net  <- Generate_BB(N = 1000, m = 15, num_seed = 100, 
                          multiple_node = 100,
                          s = 1) 
                          #shape = 1, rate = 1)
  
      net_stats <- GetStatistics(net$graph) 
  
      print(system.time(result <- JointEstimate(net$graph, 
                                                net_stats)))
      
      result
      
      
                        
      alpha_vec[i]     <- result$estimate_result$alpha
      s_vec[i]         <- result$estimate_result$shape
      alpha_optimal[i] <- result$cv_result$alpha_optimal
      r_optimal[i]     <- result$cv_result$r_optimal
  }
  print(alpha_vec)  
  mean(alpha_vec)
  print(s_vec)
  mean(s_vec)
}
