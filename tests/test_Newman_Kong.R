#testing GenerateNet, GetStatistics, Jeong, Newman_corrected and PAFit
library(PAFit)
for (prob_m in c("TRUE", "FALSE"))
  for (inc in c("TRUE","FALSE"))
    for (log in c("TRUE", "FALSE"))
      
      for (i in 1:3) {
        net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                            mode = i, shape = 100, rate = 100)
        net_stats <- GetStatistics(net$graph,deg_threshold = 1, Binning = TRUE, G = 10) 
        result_Jeong  <- Jeong(net$graph,net_stats, T_0 = 10, T_1 = 40)
        plot(result_Jeong,net_stats)
        plot(result_Jeong,net_stats, line = TRUE)
        print(result_Jeong)
        summary(result_Jeong)
        calculate_error_PA(k = result_Jeong$k,A = result_Jeong$A, mode = i)
        
        result_Jeong  <- Jeong(net$graph,net_stats, T_0 = 10, T_1 = 40, interpolate = TRUE)
        print(result_Jeong)
        summary(result_Jeong)
        plot(result_Jeong,net_stats)
        plot(result_Jeong,net_stats, line = TRUE)
        calculate_error_PA(k = result_Jeong$k,A = result_Jeong$A, mode = i)
        
        result_Newman <- Newman_corrected(net_stats)
        print(result_Newman)
        summary(result_Newman)
        plot(result_Newman,net_stats)
        plot(result_Newman,net_stats, line = TRUE)
        calculate_error_PA(k = result_Newman$k,A = result_Newman$A, mode = i)
        
        result_Newman <- Newman_corrected(net_stats, interpolate = TRUE)
        print(result_Newman)
        summary(result_Newman)
        plot(result_Newman,net_stats)
        plot(result_Newman,net_stats, line = TRUE)
        calculate_error_PA(k = result_Newman$k,A = result_Newman$A, mode = i)
}
     

