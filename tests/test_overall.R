#testing GenerateNet, GetStatistics, Jeong, Newman_corrected and PAFit
library(PAFit)
for (prob_m in c("TRUE", "FALSE"))
   for (inc in c("TRUE","FALSE"))
      for (log in c("TRUE", "FALSE"))
 
              for (i in 1:3) {
                  net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                                      mode = i, shape = 100, rate = 100)
                  for (bin in c("TRUE","FALSE"))   
                    for (deg_thresh in c(0)) {  
                        net_stats <- GetStatistics(net$graph,deg_threshold = deg_thresh, Binning = bin, G = 10) 
                        #check stats
                        if (sum(net_stats$m_t) != sum(net_stats$Sum_m_k))
                            stop("wrong at m_t and sum_m_k")
                        if (sum(abs(colSums(net_stats$m_tk) - net_stats$Sum_m_k)) != 0)
                            stop("wrong at m_tk and sum_m_k")
                        temp <- sapply(1:(net_stats$T-1),function(x) sum(net_stats$node_degree[x,] != -1))
                        if (sum(abs(rowSums(net_stats$n_tk) - (rowSums(net_stats$offset_tk) + temp))))
                            stop("wrong at node_degree, n_tk, offset_tk") 
                        if (sum(net_stats$z_j) > sum(net_stats$m_t))
                            stop("wrong at z_j")   
                        if (sum(net_stats$z_j) + sum(net_stats$offset_m_tk) - sum(net_stats$Sum_m_k))
                            stop("Wrong at offset_m_tk") 
                    }
              }
print(net_stats)
summary(net_stats)
for (mode_f_value in c("Constant_PA", "Log_linear")) {
                        result <- PAFit(net_stats, mode_f = mode_f_value,stop_cond = 10^-3)
                        print(result)
                        summary(result)
                        
                        plot(result,net_stats,plot = "A")
                        plot(result,net_stats,plot = "f")
                        plot(result,net_stats,true_f = net$fitness,plot = "true_f")
                        result <- PAFit(net_stats, 
                                        mode_f = mode_f_value,
                                        stop_cond = 10^-3, only_PA = TRUE)
                        summary(result)
                        print(result)
                        
                        result_Jeong  <- Jeong(net$graph,net_stats, T_0 = 10, T_1 = 40)
                        print(result_Jeong)
                        summary(result_Jeong)
                        plot(result_Jeong,net_stats)
                        result_Newman <- Newman_corrected(net_stats)
                        calculate_error_PA(k = result$k,A = result$A)
                        calculate_error_fitness(true = net$fitness,
                                                estimate = result$f)
      }

