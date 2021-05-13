PAFit_oneshot <- function(net,M = 100,S = 5, loop = 5,
                              G = 1000) {

result_list <- vector(mode = "list",length = loop)
for (j in 1:loop) {
  result_list[[j]] <- .estimate_multiple_phase_correct_mcPAMLE(net,M = M, S = 5, G = G,
                                                              method_name = "SI",
                                                              final_PA_type = "mean")
}


prob     <- colMeans(result_list[[1]]$final_appear)
k        <- as.integer(names(prob)) 
k_mean <- k
if (loop > 1) {
for (j in 2:loop) {
  prob <- colMeans(result_list[[j]]$final_appear)
  k    <- as.integer(names(prob)) 
  k_mean <- union(k_mean,k)
}
}

prob_matrix <- matrix(0,ncol = length(k_mean),nrow = loop)
var_matrix <-  matrix(0,ncol = length(k_mean), nrow = loop)
colnames(var_matrix) <- k_mean


colnames(prob_matrix) <- k_mean
for (j in 1:loop) {
  prob_temp <- colMeans(result_list[[j]]$final_appear)
  var_temp  <- apply(result_list[[j]]$final_appear,2,var)/dim(result_list[[j]]$final_appear)[2]
  prob_matrix[j,names(prob_temp)] <- prob_temp  
  var_matrix[j,names(var_temp)]   <- var_temp
}

input_ob <- result_list[[1]]$first_result
g  <- input_ob$g
#M  <- nrow(temp)
#print(M)
deg <- input_ob$deg
deg_bin <- input_ob$deg_bin
#M2 <- M
dist <- input_ob$dist
net_type   <- input_ob$net_type
deg_penmax <- input_ob$deg_penmax
bin_object <- input_ob$bin_object
p_estimate <- input_ob$p_estimate
Total_T    <- input_ob$Total_T
original_cumsum <- input_ob$original_cumsum
bin_dist        <- input_ob$bin_dist
bin_dist_new     <- input_ob$bin_dist_new
original_cumsum_bin <- input_ob$original_cumsum_bin
g_new <- input_ob$g_new
bin_object_new <- input_ob$bin_object_new
original_cumsum_new <- input_ob$original_cumsum_bin
deg_bin_new         <- input_ob$deg_bin


final_appear <- prob_matrix
select_appear_bin        <- rep(0,g_new)
names(select_appear_bin) <-  0:(g_new - 1) 
sd_selective_bin            <- rep(0,g_new)
names(sd_selective_bin)     <- 0:(g_new - 1) 
var_vec_select              <- rep(0,g_new)
denominator_select_bin   <- rep(0,g_new)
  
for (i in 1:(length(dist) - 1)) {
    bin_index_new <- which(bin_object_new$start <= deg[i] & deg[i] <= bin_object_new$end) # find the bin of that degree i
    #print(bin_index_new)
    if (length(bin_index_new) != 1) print("Wrong at bin index")
    denominator_select_bin[bin_index_new] <- denominator_select_bin[bin_index_new] + 
      bin_dist[i] * mean(final_appear[,deg[i] + 1])
    var_vec_select[bin_index_new] <- var_vec_select[bin_index_new] + 
      bin_dist[i]^2 * sum(var_matrix[,deg[i] + 1])/(length(var_matrix[,deg[i] + 1]))^2
}
  
sd_selective_bin       <- sqrt(var_vec_select)
  
  
selective_result_bin        <- rep(0,g_new)
names(selective_result_bin) <- 0:(g_new - 1) 
  
  
  for (i in 1:(g_new)) {
    
    if (denominator_select_bin[i] != 0) {
      selective_result_bin[i] <- original_cumsum_new[i]/denominator_select_bin[i]
    } else{
      selective_result_bin[i] <- 0
    }
  } 
  #print(selective_result_bin)
  
  sd_selective_bin[is.na(sd_selective_bin)] <- 0
  sd_selective_bin[is.nan(sd_selective_bin)] <- 0
  
  selective_result_bin[selective_result_bin == "NaN"] <- 0 
  selective_result_bin[selective_result_bin == "Inf"] <- 0 
  
  
  first_non_zero  <- which(deg_bin_new > 0)[1]
  while (selective_result_bin[first_non_zero] == 0) first_non_zero <- first_non_zero + 1
  factor_second <- selective_result_bin[first_non_zero]
  selective_result_bin <- selective_result_bin / selective_result_bin[first_non_zero]
  selective_result_bin[selective_result_bin == "NaN"] <- 0 
  selective_result_bin[selective_result_bin == "Inf"] <- 0 
  
  sd_selective_bin <- original_cumsum_new * sd_selective_bin/denominator_select_bin^2 / factor_second
  sd_selective_bin[sd_selective_bin == "NaN"] <- 0 
  sd_selective_bin[sd_selective_bin == "Inf"] <- 0 
  sd_log_sd_selective_bin <- sd_selective_bin/selective_result_bin
  sd_log_sd_selective_bin[sd_log_sd_selective_bin == "NaN"] <- 0 
  sd_log_sd_selective_bin[sd_log_sd_selective_bin == "Inf"] <- 0 
  upper <- exp(log(selective_result_bin) + 2 * sd_log_sd_selective_bin)
  lower <- exp(log(selective_result_bin) - 2 * sd_log_sd_selective_bin)
  
  
  # ignore sd, since sd does not work well here
  regress_select_bin_no_sd        <- .regress_alpha(k = bin_object_new$center_bin, 
                                                   A = selective_result_bin)
  
  
  use_var <- rep(0,g_new)
  use_var[which(sd_selective_bin == 0 & denominator_select_bin > 0)] <- 1
  #plug in the smallest possible variance
  sd_log_sd_selective_bin[use_var == 1] <- min(sd_log_sd_selective_bin[sd_log_sd_selective_bin > 0])
  
  no_zero <- bin_object_new$center_bin > 0 & selective_result_bin > 0 &
    sd_log_sd_selective_bin > 0
  if (sum(no_zero > 0) >= 2) {
    regress_select_bin <-   .regress_alpha(k = bin_object_new$center_bin, 
                                          A = selective_result_bin,
                                          sd_log_A = sd_log_sd_selective_bin)
  } else regress_select_bin <- regress_select_bin_no_sd 
  
  alpha_select_bin_no_sd <- regress_select_bin_no_sd$alpha
  alpha_select_bin <- regress_select_bin$alpha
  
  
  final_result <- list(A = selective_result_bin,
                             k = bin_object_new$center_bin,
                             alpha = alpha_select_bin,
                             alpha_no_sd = alpha_select_bin_no_sd,
                             A_upper = upper,
                             A_lower = lower,
                             regress_res = regress_select_bin,
                             regress_result_no_sd = regress_select_bin_no_sd)
final_prob <- colMeans(prob_matrix)

return(list(M = M, S= S, loop = loop,G = G,result_list = result_list,
            prob_matrix = prob_matrix, final_prob = final_prob, final_result = final_result))
}