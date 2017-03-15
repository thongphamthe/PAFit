#testing CV-related functions
library(PAFit)
prob_m <- "TRUE"
inc    <- "TRUE"
log    <-  c("TRUE")
mode_f_value <- c("Constant_PA", "Log_linear")[1]  
i <- 1
net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                    mode = i, shape = 100, rate = 100)

net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                    mode = i, shape = 100, rate = 100, specific_start = 10, 
                    mode_f = "log_normal")

net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                    mode = i, shape = 100, rate = 100, specific_start = 10, 
                    mode_f = "power_law")

net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                    mode = i, shape = 100, rate = 100, custom_PA = 1:1000,
                    specific_start = 10, 
                    mode_f = "power_law")


net_stats <- GetStatistics(net$graph,deg_threshold = 1, 
                             net_type = "directed",
                             Binning = TRUE, G = 10) 

net_stats <- GetStatistics(net$graph,deg_threshold = 0, 
                           net_type = "directed",
                           Binning = TRUE, G = 10, CompressMode = 1,
                           CompressRatio = 1, only_true_deg_matrix = TRUE)

net_stats <- GetStatistics(net$graph,deg_threshold = 0, 
                           net_type = "directed",
                           Binning = TRUE, G = 10, CompressMode = 1,
                           CompressRatio = 1)


net_stats <- GetStatistics(net$graph,deg_threshold = 0, 
                           net_type = "directed",
                           Binning = TRUE, G = 10, CompressMode = 1,
                           CompressRatio = 2)
for (n in c(20,50))
for (q in 1:3)
 for (j in 0:2)
   for (uu in 0:1) {  
      net  <- GenerateNet(N = n, m = 1,prob_m = prob_m, increase = inc, log = log,
                          mode = 1, shape = 100, rate = 100) 
      net_stats <- GetStatistics(net$graph,deg_threshold = 1, 
                                 net_type = "directed",
                                 Binning = TRUE, G = 10) 
      result <- PAFit(net_stats, r = 0.01, mode_reg_A = j, weight_PA_mode = uu,
                    stop_cond = 10^-2, q = q, normalized_f = TRUE, debug = TRUE)
      plot(result,net_stats,high_deg_A = 1)
      plot(result,net_stats)
      plot(result,net_stats,plot = "f", high_deg_f = 1)
      plot(result,net_stats,plot = "f")
   }

net_stats <- GetStatistics(net$graph,deg_threshold = 1, 
                           net_type = "directed",
                           Binning = TRUE, G = 3) 

result <- PAFit(net_stats, r = 0,  mode_reg_A = 0, s = 100,
                stop_cond = 10^-2, q = 1, normalized_f = FALSE, debug = TRUE)

result <- PAFit(net_stats, r = 0, mode_reg_A = 1, s = 100,
                stop_cond = 10^-2, q = 1, normalized_f = FALSE, debug = TRUE)

result <- PAFit(net_stats, r = 0,  mode_reg_A = 2, s = 100,
                stop_cond = 10^-2, q = 1, normalized_f = FALSE, debug = TRUE)


net_stats <- GetStatistics(net$graph,deg_threshold = 1, 
                           net_type = "directed",
                           Binning = TRUE, G = 4) 
result <- PAFit(net_stats, r = 0,  mode_reg_A = 0,
                stop_cond = 10^-2, q = 1, normalized_f = FALSE, debug = TRUE)
result <- PAFit(net_stats, r = 0, mode_reg_A = 1,
                stop_cond = 10^-2, q = 1, normalized_f = FALSE, debug = TRUE)
result <- PAFit(net_stats, r = 0,  mode_reg_A = 2,
                stop_cond = 10^-2, q = 1, normalized_f = FALSE, debug = TRUE)


print(result)
summary(result)
plot(result,net_stats)

result <- PAFit(net_stats, r = 0.01, only_PA = TRUE, 
                stop_cond = 10^-2, q = 1, debug = TRUE)
print(result)
summary(result)


result <- PAFit(net_stats, r = 0.01, only_f = TRUE, mode_f = "Log_linear",
                stop_cond = 10^-2, q = 1, normalized_f = FALSE, debug = TRUE)


result <- PAFit(net_stats, r = 0.01, mode_f = "Log_linear",
                stop_cond = 10^-2, q = 1, normalized_f = FALSE, debug = TRUE)

plot(result,net_stats,plot = "f",plot_true_degree = TRUE)

result <- PAFit(net_stats, r = 0.01, only_f = TRUE, 
                stop_cond = 10^-2, q = 1, normalized_f = FALSE, debug = TRUE)

plot(result,net_stats, plot = "f",plot_true_degree = TRUE)

plot(result,net_stats,plot = "f", high_deg_f = 2)

plot(result,net_stats,plot = "f", f_min = 0.001)
plot(result,net_stats,plot = "f", f_min = 0.001, f_max = 2)
plot(result,net_stats,plot = "f", f_max = 2)


plot(result,net_stats, true_f = net$fitness, plot = "true_f", high_deg_f = 2)

plot(result,net_stats,true_f = net$fitness, plot = "true_f" , f_min = 0.001)
plot(result,net_stats,true_f = net$fitness, plot = "true_f" , f_min = 0.001, f_max = 2)
plot(result,net_stats,true_f = net$fitness, plot = "true_f" , f_max = 2)


result <- PAFit(net_stats, r = 0.01, mode_f = "Log_linear",
                stop_cond = 10^-2, q = 1, normalized_f = FALSE, debug = TRUE)
plot(result,net_stats,line = TRUE)




print(result)
summary(result)


result <- PAFit(net_stats, r = 0.01, only_f = TRUE, 
                stop_cond = 10^-2, q = 1, normalized_f = TRUE, debug = TRUE)


net_stats <- GetStatistics(net$graph,deg_threshold = 0, 
                           net_type = "directed",
                           Binning = TRUE, G = 1, CompressMode = 0, CompressRatio = 0.5)

result <- PAFit(net_stats, r = 0.01, 
                stop_cond = 10^-2, q = 1, normalized_f = TRUE, debug = TRUE, start_mode_A = "Random")

result <- PAFit(net_stats, r = 0.01, only_f = TRUE,
                stop_cond = 10^-2, q = 1, normalized_f = TRUE, debug = TRUE, true_A = 1:10000)


net_stats <- GetStatistics(net$graph,deg_threshold = 0, 
                           net_type = "directed",
                           Binning = TRUE, G = 10, CompressMode = 0, 
                           CompressRatio = 2)
result <- PAFit(net_stats, r = 0.01, 
                stop_cond = 10^-2, q = 1, normalized_f = TRUE, debug = TRUE, start_mode_A = "Random")

result <- PAFit(net_stats, r = 0.01, only_f = TRUE,
                stop_cond = 10^-2, q = 1, normalized_f = TRUE, debug = TRUE, true_A = 1:10000)





