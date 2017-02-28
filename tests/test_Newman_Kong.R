#testing GenerateNet, GetStatistics, Jeong, Newman_corrected and PAFit

library(PAFit)
for (prob_m in c("TRUE", "FALSE"))
  for (inc in c("TRUE","FALSE"))
    for (log in c("TRUE", "FALSE"))
      
      for (i in 1:3) {
        net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                            mode = i, shape = 100, rate = 100)
        net_stats <- GetStatistics(net$graph,deg_threshold = 1, Binning = TRUE, G = 10) 
        result_Jeong  <- Jeong(net$graph , net_stats, T_0_start = 0, T_0_end = 20, T_1_start = 30)
        plot(result_Jeong,net_stats)
        plot(result_Jeong,net_stats, line = TRUE)
        print(result_Jeong)
        summary(result_Jeong)
       
        result_Jeong  <- Jeong(net$graph,net_stats, T_0_start = 0, T_0_end = 20, T_1_start = 30)
        print(result_Jeong)
        summary(result_Jeong)
        plot(result_Jeong,net_stats)
        plot(result_Jeong,net_stats, line = TRUE)
       
        result_Newman <- Newman(net_stats)
        print(result_Newman)
        summary(result_Newman)
        plot(result_Newman,net_stats)
        plot(result_Newman,net_stats, line = TRUE)
       result_Newman <- Newman(net_stats, interpolate = TRUE)
        print(result_Newman)
        summary(result_Newman)
        plot(result_Newman,net_stats)
        plot(result_Newman,net_stats, line = TRUE)
       
}
     
