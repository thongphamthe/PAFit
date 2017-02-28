#testing CV-related functions
library(PAFit)
prob_m <- "TRUE"
inc    <- "TRUE"
log    <-  c("TRUE")
mode_f_value <- c("Constant_PA", "Log_linear")[1]  
i      <- 1
net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                    mode = i, shape = 100, rate = 100,alpha = 0.5)



net_stats <- GetStatistics(net$graph,deg_threshold = 1, 
                           net_type = "directed",
                           Binning = TRUE, G = 10) 

result <- PAFit(net_stats, mode_f = "Log_linear", 
                true_f = rep(1,length(net_stats$f_position)),
                s = 100, debug = TRUE)

result <- PAFit(net_stats, only_f = TRUE,
                true_A = rep(1,net_stats$G),
                s = 100, debug = TRUE)



print(result)
