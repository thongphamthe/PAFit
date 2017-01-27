#testing CV-related functions
library(PAFit)
for (prob_m in c("TRUE", "FALSE"))
  for (inc in c("TRUE","FALSE"))
    for (log in c("TRUE", "FALSE"))
      for (mode_f_value in c("Constant_PA", "Log_linear"))  
        for (i in 1:3) {
          net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                              mode = i, shape = 100, rate = 100)
          for (bin in c("TRUE","FALSE"))   
            for (net_type in c("directed","undirected")) {  
              net_stats <- GetStatistics(net$graph,deg_threshold = 1, 
                                         net_type = net_type,
                                         Binning = bin, G = 10) 
              
              # create CV_data
              data_cv <- CreateDataCV(net$graph, p = 0.75, G = 50, 
                           net_type = net_type,deg_thresh = 0)
              # perform CV  
              cv_result <- performCV(data_cv,stop_cond = 10^-3)
            }
}
