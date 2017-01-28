#testing CV-related functions
library(PAFit)
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
for (q in 1:3)
 for (j in 0:2)
   for (uu in 0:1)   
      result <- PAFit(net_stats, r = 0.01, mode_reg_A = j, weight_PA_mode = uu,
                    stop_cond = 10^-5, q = q, normalized_f = TRUE)

result <- PAFit(net_stats, r = 0, estimate_shape = TRUE, 
                stop_cond = 10^-5, q = 1, normalized_f = FALSE)
print(result)
summary(result)

result <- PAFit(net_stats, r = 0.01, only_PA = TRUE, 
                stop_cond = 10^-5, q = 1)
print(result)
summary(result)


result <- PAFit(net_stats, r = 0.01, only_f = TRUE, mode_f = "Log_linear",
                stop_cond = 10^-5, q = 1, normalized_f = FALSE)


result <- PAFit(net_stats, r = 0.01, mode_f = "Log_linear",
                stop_cond = 10^-5, q = 1, normalized_f = FALSE)



result <- PAFit(net_stats, r = 0.01, only_f = TRUE, 
                stop_cond = 10^-5, q = 1, normalized_f = FALSE)


print(result)
summary(result)


result <- PAFit(net_stats, r = 0.01, only_f = TRUE, estimate_shape = TRUE, 
                stop_cond = 10^-5, q = 1, normalized_f = TRUE)

