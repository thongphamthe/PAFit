#testing CV-related functions
library(PAFit)
prob_m <- "TRUE"
inc    <- "TRUE"
log    <-  c("TRUE")
mode_f_value <- c("Constant_PA", "Log_linear")[1]  
i      <- 1
net  <- generate_net(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                    mode = i, s = 100,alpha = 0.5)



net_stats <- get_statistics(net,deg_threshold = 1, 
                           binning = TRUE, g = 10) 



