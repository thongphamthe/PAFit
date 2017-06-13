#testing CV-related functions

## code here that contains the tests
  

rm(list = ls())
library(PAFit)
set.seed(1)
prob_m <- "FALSE"
inc    <- "FALSE"
log    <-  c("FALSE")
i      <- 1

net  <- generate_net(N = 1000, m = 5,prob_m = prob_m, 
                    increase = inc, log = log, multiple_node = 100, num_seed = 100,
                    mode = i, shape = 10, rate = 10,alpha = 1)

net_stats <- get_statistics(net,deg_threshold = 1, 
                           binning = TRUE, g = 50) 

result <- joint_estimate(net, net_stats)
plot(result,net_stats)




