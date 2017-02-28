#testing CV-related functions

library(PAFit)
prob_m <- "TRUE"
inc    <- "TRUE"
log    <-  c("TRUE")
mode_f_value <- c("Constant_PA", "Log_linear")[1]  
i      <- 1
net  <- GenerateNet(N = 500, m = 50,prob_m = prob_m, increase = inc, log = log,
                    mode = i, shape = 100, rate = 100,alpha = 0.5)



net_stats <- GetStatistics(net$graph,deg_threshold = 1, 
                           net_type = "directed",
                           Binning = TRUE, G = 10) 

result <- PAFit(net_stats, only_f = TRUE, start_f = rep(1,length(net_stats$f_position)),
                true_A = 1:net_stats$G,s = 100, debug = TRUE, stop_cond = 10^-8)


#print(result)
