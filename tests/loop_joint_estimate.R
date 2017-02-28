if (FALSE) {
  #setwd("tests")
  rm(list = ls())
  library(PAFit)
  set.seed(1)
  prob_m <- "FALSE"
  inc    <- "FALSE"
  log    <-  c("FALSE")
  
  M <- 20
  
  alpha_vec <- rep(0,M)
  s_vec     <- rep(0,M)
  alpha_optimal <- rep(0,M)
  
  for (i in 1:M) {
  
      net  <- GenerateNet(N = 2000, m = 50,prob_m = prob_m, num_seed = 400, multiple_node = 200,
                          increase = inc, log = log,
                          mode = 1, shape = 2, rate = 2,alpha = 1, noNewNodeStep = 5, m_noNewNodeStep = 2000)
  
      net_stats <- GetStatistics(net$graph,deg_threshold = 0, 
                                net_type = "directed",
                                 Binning = TRUE, G = 50) 
  
      result <- JointEstimate(raw_net = net$graph, net_stat = net_stats, mode_reg_A = 0)
      
      alpha_vec[i] <- result$estimate_result$alpha
      s_vec[i]     <- result$estimate_result$shape
      alpha_optimal[i] <- result$cv_result$alpha_optimal
  }
  print(alpha_vec)  
  print(s_vec)
}