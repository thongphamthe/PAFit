#testing codes
result <- NULL
data   <- NULL
stats  <- NULL
library(PAFit)
#Repeat the test 5 times 
for (count in NULL)
for (prob_m in c("TRUE", "FALSE")[1])
   for (inc in c("TRUE","FALSE")[1])
      for (log in c("TRUE", "FALSE")[1])             
              for (i in 1:1) {
                  data  <- GenerateNet(N = 10, m = 10,prob_m = prob_m, increase = inc, log = log,
                                      mode = i, shape = 2, rate = 2)
                  for (bin in c("TRUE","FALSE")[1])   
                    for (deg_thresh in c(0,2)[1]) {  
                        stats <- GetStatistics(data$graph,deg_threshold = deg_thresh, Binning = bin, G = 10) 
                        #check stats
                        if (sum(stats$m_t) != sum(stats$Sum_m_k))
                            print("wrong at m_t and sum_m_k")
                        if (sum(abs(colSums(stats$m_tk) - stats$Sum_m_k)) != 0)
                            print("wrong at m_tk and sum_m_k")
                        temp <- sapply(1:(stats$T-1),function(x) sum(stats$node_degree[x,] != -1))
                        if (sum(abs(rowSums(stats$n_tk) - (rowSums(stats$offset_tk) + temp))))
                            print("wrong at node_degree, n_tk, offset_tk") 
                        if (sum(stats$z_j) > sum(stats$m_t))
                            print("wrong at z_j")   
                        if (sum(stats$z_j) + sum(stats$offset_m_tk) - sum(stats$Sum_m_k))
                            print("Wrong at offset_m_tk")  
                        result <- PAFit(stats, mode_f ="Constant_PA",stop_cond = 10^-3)
                        plot(result,stats,plot = "A")
                        #plot(result,stats,plot = "f")
                        #plot(result,data = stats,true_f = data$fitness,plot = "true_f")
                  }
    }
