#testing generate_net, get_statistics, Jeong, Newman_corrected and PAFit

library(PAFit)
for (prob_m in c("TRUE", "FALSE"))
   for (inc in c("TRUE","FALSE"))
      for (log in c("TRUE", "FALSE"))
 
              for (i in 1:3) {
                  net  <- generate_net(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                                      mode = i, s = 10)
                  for (bin in c("TRUE","FALSE"))
                    
                    for (deg_thresh in c(0)) {  
                        net_stats <- get_statistics(net,deg_threshold = deg_thresh, binning = bin, g = 10) 
                        #check stats
                        if (sum(net_stats$m_t) != sum(net_stats$sum_m_k))
                            stop("wrong at m_t and sum_m_k")
                        if (sum(abs(colSums(net_stats$m_tk) - net_stats$sum_m_k)) != 0)
                            stop("wrong at m_tk and sum_m_k")
                        temp <- sapply(1:(net_stats$T-1),function(x) sum(net_stats$node_degree[x,] != -1))
                        if (sum(abs(rowSums(net_stats$n_tk) - (rowSums(net_stats$offset_tk) + temp))))
                            stop("wrong at node_degree, n_tk, offset_tk") 
                        if (sum(net_stats$z_j) > sum(net_stats$m_t))
                            stop("wrong at z_j")   
                        if (sum(net_stats$z_j) + sum(net_stats$offset_m_tk) - sum(net_stats$sum_m_k))
                            stop("Wrong at offset_m_tk") 
                        net_stats <- get_statistics(net,deg_threshold = deg_thresh, binning = bin, g = 10,only_PA = TRUE)
                        net_stats <- get_statistics(net,deg_threshold = deg_thresh, binning = bin, g = 10,only_true_deg_matrix = TRUE)
                    }
              }
net_stats <- get_statistics(net,deg_threshold = deg_thresh, binning = TRUE, g = 10) 
print(net_stats)
summary(net_stats)
for (mode_f_value in c("Constant_PA", "Log_linear")) {
                        result_Jeong  <- Jeong(net,net_stats, T_0_start = 0, T_0_end = 20, T_1_start = 30 , T_1_end = 40)
                        result_Jeong  <- Jeong(net,net_stats, T_0_start = 0, T_0_end = 20, T_1_start = 30 , T_1_end = 40, interpolate = TRUE)
                        print(result_Jeong)
                        summary(result_Jeong)
                        plot(result_Jeong,net_stats)
                        plot(result_Jeong,net_stats, line = TRUE)
                        plot(result_Jeong,net_stats, high_deg = 5)
                        result_Newman <- Newman(net, net_stats)
                        result_Newman <- Newman(net, net_stats, interpolate = TRUE)
                      
}
net_stats <- get_statistics(net,deg_threshold = deg_thresh, binning = FALSE, g = 10) 
print(net_stats)
summary(net_stats)
for (mode_f_value in c("Constant_PA", "Log_linear")) {

  result_Jeong  <- Jeong(net,net_stats, T_0_start = 0, T_0_end = 20, T_1_start = 30 , T_1_end = 40)
  result_Jeong  <- Jeong(net,net_stats, T_0_start = 0, T_0_end = 20, T_1_start = 30 , T_1_end = 40, interpolate = TRUE)
  print(result_Jeong)
  summary(result_Jeong)
  plot(result_Jeong,net_stats)
  plot(result_Jeong,net_stats, line = TRUE)
  plot(result_Jeong,net_stats, high_deg = 5)
  result_Newman <- Newman(net, net_stats)
  result_Newman <- Newman(net, net_stats, interpolate = TRUE)
}

