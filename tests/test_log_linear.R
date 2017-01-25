#testing codes
result <- NULL
data   <- NULL
stats  <- NULL
library(PAFit)

data  <- GenerateNet(N = 100, m = 3, num_seed = 20, shape = 5, alpha = 0.7)
stats <- GetStatistics(data$graph, Binning = FALSE) 
result <- PAFit(stats,mode_f = "Log_linear", s = 5, debug = TRUE, stop_cond = 10^-4)

