#testing codes
result <- NULL
data   <- NULL
stats  <- NULL
library(PAFit)

net       <- GenerateNet(N = 100, m = 1, num_seed = 20, shape = 5, alpha = 0.7)
net_stats <- GetStatistics(net$graph, Binning = FALSE) 
result    <- PAFit(net_stats,mode_f = "Log_linear", s = 5, debug = TRUE, stop_cond = 10^-4)

