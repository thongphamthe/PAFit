#testing CV-related functions

## code here that contains the tests


rm(list = ls())
library(PAFit)


net_object  <- generate_net(N = 20, m = 10,
                            s = 50,alpha = 0.5)
net_stat <- get_statistics(net_object) 
result <- joint_estimate(net_object,net_stat, stop_cond = 10^-3)

simulated_data <- generate_simulated_data_from_estimated_model(net_object,net_stat,result, M = 3)

plot_contribution(simulated_data, result, which_plot = "PA")

plot_contribution(simulated_data, result, which_plot = "fit")

