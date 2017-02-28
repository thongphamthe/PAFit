if (FALSE) {
  #setwd("tests")
  rm(list = ls())
  library(PAFit)
  set.seed(1)
  prob_m <- "FALSE"
  inc    <- "FALSE"
  log    <-  c("FALSE")
  

    
    net  <- GenerateNet(N = 1000, m = 50,prob_m = prob_m, num_seed = 100, multiple_node = 100,
                        increase = inc, log = log,
                        mode = 1, shape = 100, rate = 100,alpha = 0.5)
    
    net_stats <- GetStatistics(net$graph,deg_threshold = 0, 
                               net_type = "directed",
                               Binning = TRUE, G = 50) 
    
    result <- PAFit_new(net_stats,r = 1, s = 100, mode_f = "Log_linear")
    
    plot(result,net_stats)
    
    plot(result,net_stats,true_f = net$fitness, plot = "true_f", high_deg_f = 1)
    
    
    T <- max(net_stats$appear_time)
    var_array <- rep(0,T)
    plot("", xlim = c(min(result$f),max(result$f)), ylim = c(1,T), xlab =  "Fit", ylab = "Time")
    for (i in 1:T) {
      f <- result$f[net_stats$appear_time == i]  
      var_array[i] <- var(f)
       points(f,rep(i,length(f)))  
    }
    
    
}