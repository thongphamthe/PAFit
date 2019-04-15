#testing CV-related functions

## code here that contains the tests
  

rm(list = ls())
library(PAFit)
# alpha = 0.5, s = 10, ii = 220
# for CRAN. In developing, set ii from 1 to 1000
for (ii in 1) {
    set.seed(2)
    print(ii)
    prob_m <- "FALSE"
    inc    <- "FALSE"
    log    <-  c("FALSE")
    i      <- 1

    net  <- generate_net(N = 100, m = 15,prob_m = prob_m, 
                         increase = inc, log = log, multiple_node = 20, num_seed = 20,
                         mode = i, s = 10,alpha = 0.5)

    net_stats <- get_statistics(net,deg_threshold = 5, 
                                binning = TRUE, g = 50) 

    result_A <- only_A_estimate(net, net_stats, stop_cond = 10^-5)
    
    result <- joint_estimate(net, net_stats, stop_cond = 10^-5)
    plot(result,net_stats)
}




