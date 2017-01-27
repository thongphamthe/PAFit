#testing CV-related functions
library(PAFit)
prob_m <- "TRUE"
inc    <- "TRUE"
log <-  c("TRUE")
mode_f_value <- c("Constant_PA", "Log_linear")[1]  
i <- 1
net  <- GenerateNet(N = 50, m = 10,prob_m = prob_m, increase = inc, log = log,
                   mode = i, shape = 100, rate = 100)
for (net_type in c("directed","undirected")) {  
    net_stats <- GetStatistics(net$graph,deg_threshold = 1, 
                               net_type = net_type,
                               Binning = TRUE, G = 10) 
              # create CV_data
    data_cv <- CreateDataCV(net$graph, p = 0.75, G = 50, 
                            net_type = net_type,deg_thresh = 0)
}
#print(data_cv)
#summary(data_cv)

   # perform CV  
cv_result <- performCV(data_cv,  r = 10^c(-2,0),s = 10^c(-1,1),
                       stop_cond = 10^-3,only_PAFit = FALSE)
#print(cv_result)
#summary(cv_result)


cv_result <- performCV(data_cv,  r = 10^c(-2,0),s = 10^c(-1,1),
                       stop_cond = 10^-3,only_PAFit = TRUE)


cv_result <- performCV(data_cv, r = 10^c(-2,0),s = 10^c(-1,1),
                       stop_cond = 10^-3, only_linear = TRUE)

