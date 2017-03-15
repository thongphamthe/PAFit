#testing CV-related functions

## code here that contains the tests
  

rm(list = ls())
library(PAFit)
set.seed(1)
prob_m <- "FALSE"
inc    <- "FALSE"
log    <-  c("FALSE")


i      <- 1

net  <- GenerateNet(N = 20, m = 50,prob_m = prob_m, 
                    increase = inc, log = log,
                    mode = i, shape = 10, rate = 10,alpha = 0.5)

net_stats <- GetStatistics(net$graph,deg_threshold = 0, 
                           net_type = "directed",
                           Binning = TRUE, G = 50) 

result <- JointEstimate(raw_net = net$graph, 
                        net_stat = net_stats, 
                        stop.cond = 10^-3,
                         print.out = TRUE)


#print(result)
