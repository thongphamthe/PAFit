.estimate_multiple_phase_correct_mcPAMLE <- function(net, 
                                    M = 100, 
                                    S = 2,
                                    net_type    = "directed",
                                    true_time   = NULL,
                                    final_PA_type = "mean",
                                    method_name = c("SI"),
                                    G = 1000) {
  if (S < 1)
    stop("S must be at least 1.")
  
  result_1          <- .estimate_new_correct_mcPAMLE(net,M = M, net_type = net_type, 
                                                    final_PA_type = final_PA_type,
                                                    true_time = true_time, G = G);
  appear_time_list  <- vector(mode = "list", length = S)
  final_appear_list <- vector(mode = "list", length = S)
  estimated_A       <- vector(mode = "list",length = S)
  estimated_k       <- vector(mode = "list", length = S)
  est_A_upper       <- vector(mode = "list", length = S)
  est_A_lower       <- vector(mode = "list", length = S)
  detailed_result   <- vector(mode = "list", length = S)
  estimated_alpha   <- rep(0,length = S)
  
  appear_time_list[[1]]     <- result_1$appear_second
  final_appear_list[[1]]    <- result_1$final_appear#/result_1$final_appear[2]
  
  {
    estimated_A[[1]]     <- result_1$selective_result_with_bin$A
    est_A_upper[[1]]     <- result_1$selective_result_with_bin$A_upper
    est_A_lower[[1]]     <- result_1$selective_result_with_bin$A_lower
    estimated_k[[1]]     <- result_1$selective_result_with_bin$k
    estimated_alpha[1] <- result_1$selective_result_with_bin$alpha
    detailed_result[[1]] <- result_1$selective_result_with_bin
  }
  
  if (S > 1) {
    for (jj in 2:S) {
      if (jj == 2) {
       if ("SI" == method_name) {
          input_A <- result_1$selective_result_with_bin$A   
        }
      } else { 
        if ("SI" == method_name) {
          input_A <- result$selective_result_with_bin$A   
        }
      }
      
      result <- estimate_single_phase_correct_mcPAMLE(result_1,input_A,
                                                  type = method_name, final_PA_type = final_PA_type)
      
      {
        estimated_A[[jj]]     <- result$selective_result_with_bin$A
        est_A_upper[[jj]]     <- result$selective_result_with_bin$A_upper
        est_A_lower[[jj]]     <- result$selective_result_with_bin$A_lower
        estimated_k[[jj]]     <- result$selective_result_with_bin$k
        estimated_alpha[jj]   <- result$selective_result_with_bin$alpha
        detailed_result[[jj]] <- result$selective_result_with_bin
        appear_time_list[[jj]]     <- result$appear_second
        final_appear_list[[jj]]    <- result$final_appear#/result$final_appear[2]
      }
    }
  }
  
  # two way to calculate the average of A:
  # first: pool together S estimations of A
  temp      <- matrix(unlist(estimated_A),nrow = S, byrow = TRUE)
  A_average <- colMeans(temp)
  sd_A      <- apply(temp,MARGIN = 2,sd)/sqrt(S)
  sd_log_A  <- sd_A/A_average
  
  regress_temp <- .regress_alpha(k = estimated_k[[1]], 
                                A = A_average,
                                sd_log_A = sd_log_A) 
  
  first_way_average <- list(k = estimated_k[[1]],
                            A = A_average,
                            sd_A = sd_A,
                            sd_log_A = sd_log_A,
                            alpha = regress_temp$alpha,
                            regree_res = regress_temp)
  
  # second: pool together the statistics
    if ("SI" == method_name)  {
      temp <- final_appear_list[[1]]
      if (S > 1) {
        for (jjj in 1:S)
          temp <- rbind(temp,final_appear_list[[jjj]])
      }
    }
  
  input_ob <- result_1
  g  <- input_ob$g
  M  <- nrow(temp)
  #print(M)
  deg <- input_ob$deg
  deg_bin <- input_ob$deg_bin
  M2 <- M
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
  
 if ("SI" == method_name) {
    final_appear <- temp 
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
                                       var(bin_dist[i] * final_appear[,deg[i] + 1])/(length(final_appear[,deg[i] + 1]))
    }
    #print(select_mat_temp)
  
    sd_selective_bin       <- sqrt(var_vec_select)
    
   #print(sd_selective_bin)
    #print(use_var)
    
    #print(denominator_select_bin)
    #print(sd_selective_bin)
    
    
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
    }
    else regress_select_bin <- regress_select_bin_no_sd 
    
    alpha_select_bin_no_sd <- regress_select_bin_no_sd$alpha
    alpha_select_bin <- regress_select_bin$alpha
    
    
    second_way_average <- list(A = selective_result_bin,
                               k = bin_object_new$center_bin,
                               alpha = alpha_select_bin,
                               alpha_no_sd = alpha_select_bin_no_sd,
                               A_upper = upper,
                               A_lower = lower,
                               regress_res = regress_select_bin,
                               egress_result_no_sd = regress_select_bin_no_sd)
    
  }
  
  #final_appear <- rbind(final_appear,rep(1,dim(final_appear)[2]))
  return(list(first_result = result_1,
              estimated_A = estimated_A,
              estimated_k = estimated_k,
              estimated_alpha = estimated_alpha,
              detailed_result = detailed_result,
              method_name = method_name,
              appear_time_list = appear_time_list,
              final_appear_list = final_appear_list, 
              M = M,
              S = S,
              net_type = net_type,
              true_time = true_time,
              final_appear = final_appear,
              first_way_average = first_way_average,
              second_way_average = second_way_average))
}
