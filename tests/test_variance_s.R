if (FALSE) {
  #setwd("tests")
  rm(list = ls())
  library(PAFit)
  set.seed(1)
  prob_m <- "FALSE"
  inc    <- "FALSE"
  log    <-  c("FALSE")
  

    
    net  <- generate_net(N = 1000, m = 50,prob_m = prob_m, num_seed = 100, multiple_node = 100,
                        increase = inc, log = log,
                        mode = 1, s = 100,alpha = 0.5)
    
    net_stats <- get_statistics(net,deg_threshold = 0, 
                               binning = TRUE, g = 50) 
    
  
    
    T <- max(net_stats$appear_time)
    var_array <- rep(0,T)
    plot("", xlim = c(min(result$f),max(result$f)), ylim = c(1,T), xlab =  "Fit", ylab = "Time")
    for (i in 1:T) {
      f <- result$f[net_stats$appear_time == i]  
      var_array[i] <- var(f)
       points(f,rep(i,length(f)))  
    }
    
    
}