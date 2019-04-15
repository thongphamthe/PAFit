#testing CV-related functions


rm(list = ls())
library(PAFit)
# for CRAN. In developing, set ii from 1 to 1000
for (ii in 1) {
    set.seed(2)
    print(ii)
    prob_m <- "FALSE"
    inc    <- "FALSE"
    log    <-  c("FALSE")
    i      <- 1

    net  <- generate_net(N = 1000, m = 5,prob_m = prob_m, 
                         increase = inc, log = log, multiple_node = 100, num_seed = 100,
                         mode = i, s = 0,alpha = 1)

    net_stats <- get_statistics(net,deg_threshold = 1, 
                                binning = TRUE, g = 50) 
    result <- only_A_estimate(net, net_stats, stop_cond = 10^-5)
}
