#testing CV-related functions
#library(PAFit)
prob_m <- "TRUE"
inc    <- "TRUE"
log    <-  c("TRUE")
mode_f_value <- c("Constant_PA", "Log_linear")[1]  
i <- 1
net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                    mode = i, shape = 100, rate = 100)



net_stats <- GetStatistics(net$graph,deg_threshold = 1, 
                           net_type = "directed",
                           Binning = TRUE, G = 10) 

result <- PAFit(net_stats,mode_f = "Log_linear",
                auto_stop =  TRUE, s = 0.1,
                stop_cond = 10^-5,normalized_f = FALSE, debug = TRUE)
