if (FALSE) {
rm(list = ls())
set.seed(1)
library("PAFit")
# size of initial network = 100
# number of new nodes at each time-step = 100
# Ak = k; inverse variance of the distribution of node fitnesse = 5

#alpha <- rep(0,100)
s <- rep(0,50)
alpha <- rep(0,50)
for (j in 1:50) {
print(paste0("j:",j))
net        <- generate_net(N        = 1000 , m             = 50 , 
                          num_seed = 200  , multiple_node = 50, alpha = -1,
                          s        = 5)
net_stats  <- get_statistics(net,deg_threshold = 0)
result_big <- joint_estimate(net,net_stats)
s[j] <- result_big$estimate_result$shape
alpha[j] <- result_big$estimate_result$alpha
#data_cv <- .CreateDataCV(net, g = 50, deg_thresh = 0,
#                          p = 0.75)

#cv_result <- PAFit(data_cv$stats,s = 5,mode_f = "Log_linear")
#alpha[j] <- cv_result$alpha
}
}
# Joint estimation of attachment function Ak and node fitness
#result     <- joint_estimate(net, net_stats)


#data_cv <- .CreateDataCV(net, g = 50, deg_thresh = 0,
#                          p = 0.75)

#cv_result <- PAFit(data_cv$stats, r = 0.1, s = 5)

# name_z_j <- intersect(names(data_cv$stats$z_j), as.character(as.numeric(data_cv$stats$node_before_final)))
# z_j_before_final <- data_cv$stats$z_j[name_z_j]
# no_edge <- names(z_j_before_final[z_j_before_final == 0])
# cv_result$f[no_edge]
# cv_result$var_f[no_edge]
# plot(cv_result,data_cv$stats)
# cv_result




#result_only_F <- only_F_estimate(net,net_stats)

# plot(result_big,net_stats)
# plot(result_big, net_stats, true = net$fitness, plot = "true_f")
# 
