if (FALSE) {
set.seed(1)
library("PAFit")
# size of initial network = 100
# number of new nodes at each time-step = 100
# Ak = k; inverse variance of the distribution of node fitnesse = 5
net        <- generate_BB(N        = 1000 , m             = 50 , 
                          num_seed = 100  , multiple_node = 100,
                          s        = 5)
net_stats  <- get_statistics(net,deg_threshold = 5)

# Joint estimation of attachment function Ak and node fitness
result     <- joint_estimate(net, net_stats)


#data_cv <- .CreateDataCV(net, g = net_stats$g, deg_thresh = 0,
#                          p = 0.75)

#cv_result <- PAFit(data_cv$stats, mode_f = "Log_linear", r = 0.01, s = 0.1)


summary(result)

# plot the estimated attachment function
plot(result , net_stats)

# true function
true_A     <- pmax(result$estimate_result$center_k,1)
lines(result$estimate_result$center_k, true_A, col = "red") # true line
legend("topleft" , legend = "True function" , col = "red" , lty = 1 , bty = "n")

# plot the estimated node fitnesses and true node fitnesses
plot(result, net_stats, true = net$fitness, plot = "true_f")

name_z_j <- intersect(as.character(as.numeric(net_stats$node_before_final)),names(net_stats$z_j))
z_j_before_final <- net_stats$z_j[name_z_j]
no_edge <- names(z_j_before_final[z_j_before_final == 0])
result$estimate_result$f[no_edge]
result$estimate_result$var_f[no_edge]



test_result <- PAFit(net_stats, r = 1, s = 5)

z_j_before_final <- net_stats$z_j[as.character(as.numeric(net_stats$node_before_final))]
no_edge <- names(z_j_before_final[z_j_before_final == 0])
test_result$f[no_edge]
test_result$var_f[no_edge]
plot(test_result$objective_value)
}
