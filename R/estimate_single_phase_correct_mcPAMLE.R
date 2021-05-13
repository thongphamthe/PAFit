
# input: some estimation result (e.g., from estimate_new)
# output: do one time loop
# to save time, should do one time loop for both mcPAFIT-MLE and mcPAFIT-SI

# input_ob: returned object of estimate_new
estimate_single_phase_correct_mcPAMLE <- function(input_ob, first_result = input_ob$second_result_with_bin$A,
                                              type = c("MLE","SI"), final_PA_type = "mean") {
  
  g  <- input_ob$g
  M  <- input_ob$M
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
  degree_exist        <- input_ob$degree_exist
  
  
  
  is_directed <- 1
  if (net_type == "directed") {
    is_directed <- 1
  } else if (net_type == "undirected") {
    is_directed <- 0
  } else {
    stop("Wrong")
  }
  
  
  second_result    <- rep(0,length(first_result))
  
  # create a custom_PA vector by filling in the zeros
  custom_PA <- rep(1,deg_penmax + 1) # PA value from 0 to deg_penmax
  for (i in 1:(deg_penmax + 1)) {
    bin_index <- which(bin_object$start <= i - 1 & bin_object$end >= i - 1)
    if (first_result[bin_index] > 0) { custom_PA[i] <- first_result[bin_index];
    } else { 
      # if bin_index is at either end, then increase/decrease the bin index until the first value that is non-zero  
      if (bin_index == 1) {
        while (first_result[bin_index] == 0) bin_index <- bin_index + 1;
        custom_PA[i] <- first_result[bin_index];
      } else if (bin_index == g) {
        while (first_result[bin_index] == 0) bin_index <- bin_index - 1;
        custom_PA[i] <- first_result[bin_index];
      } else { 
        #interpolate on log-scale
        before <- bin_index - 1;
        # decrease the before index until the first non-zero value is reached
        while (first_result[before] == 0) {before <- before - 1; if (before <= 0) break;}
        
        after  <- bin_index + 1;
        # increase the after index until the first non-zero value is reached
        if (after <= length(first_result)) {
          while (first_result[after] == 0) {after <- after + 1; if (after > length(first_result)) break;}
        }
        if (before >0 & after <= length(first_result)) {
          expo <- log(first_result[before]/first_result[after])/log((deg_bin[before] + 1)/(deg_bin[after] + 1))
          a    <- first_result[before] / (deg_bin[before] + 1)^expo
          #if ("LogLinear" == which_interpolate) {
          #  custom_PA[i] <- a * (deg_bin[i] + 1)^expo
          #} else if ("GeoMean" == which_interpolate) {
          custom_PA[i]  <- sqrt(first_result[before] * first_result[after])
          #}
        } else {
          if (before > 0) {custom_PA[i] <- first_result[before]
          } else if (after <= length(first_result)) {
            custom_PA[i] <- first_result[after]
          }  
        }
      }
    }
  }  
  
  
  #first_non_zero <- 1
  first_non_zero  <- which(deg_bin == 1)
  while (custom_PA[first_non_zero] == 0) first_non_zero <- first_non_zero + 1;
  custom_PA <- custom_PA / custom_PA[first_non_zero]
  
  
  if ("mean" == final_PA_type) {
    final_PA_val <- mean(custom_PA[custom_PA > 0])  
  } else if ("max" == final_PA_type) {
    final_PA_val <- max(custom_PA[custom_PA > 0])  
  }
  appear_second <- matrix(0,nrow = M2, ncol = g)
  colnames(appear_second) <- 0:(g-1)
  final_appear <- matrix(0,nrow = M2, ncol = g)
  colnames(final_appear) <- 0:(g-1)
  
  mode_PA    = 3;
  bin_vector = bin_object$bin
  alpha = 1;
  beta =  2;
  #print(dim(final_appear))
  #print(g)
  .generate_net_C_with_count_multi_corrected(appear_second,
                                   final_appear,
                                   Total_T, mode_PA, alpha,beta, custom_PA, p_estimate,
                                   bin_vector , g, 
                                   deg_penmax,final_PA_val,is_directed,M,
                                   degree_exist)
  
  #print(dim(final_appear))
  ## prepare the variables
  alpha_select_bin_no_sd     <- NULL
  alpha_select_bin           <- NULL
  selective_result_with_bin  <- NULL
  regress_res_bin            <- NULL
  alpha_second_bin           <- NULL
  second_result_with_bin     <- NULL
  sd_selective_bin           <- NULL
  sd_result_second_bin       <- NULL
  
  ### mcPAFIT-MLE ####
  if ("MLE" == type) {
    {
      second_result_bin        <- rep(0,g_new)
      names(second_result_bin) <- 0:(g_new-1)
      

      
      
      denominator_bin          <- rep(0,g_new)
      names(denominator_bin)   <-  0:(g_new - 1) 
      
      var_temp                 <- rep(0,g_new)
      names(var_temp)          <- 0:(g_new-1)
      
      
      for (i in 1:(length(dist) - 1)) {
        bin_index_new <- which(bin_object_new$start <= deg[i] & deg[i] <= bin_object_new$end) # find the bin of that degree i
        
        appear_vec_sec <- appear_second[,deg[i] + 1] 
        non_negative   <- appear_vec_sec >= 0 
        if (sum(appear_vec_sec >0) > 0) {
          use_vec                 <- appear_vec_sec[non_negative]
          var_temp[bin_index_new] <- var_temp[bin_index_new] + var(use_vec*bin_dist[i])/(length(use_vec))
          denominator_bin[bin_index_new] <- denominator_bin[bin_index_new] + bin_dist[i] * mean(use_vec)   
        }
      }
      sd_second_bin                <- sqrt(var_temp)
      names(sd_second_bin)         <- 0:(g_new-1)
      sd_second_bin[is.na(sd_second_bin)] <- 0
      sd_second_bin[is.nan(sd_second_bin)] <- 0
      
      for (i in 1:(g_new)) {
        if (is.na(sd_second_bin[i])) sd_second_bin[i] <- 0
        if (denominator_bin[i] != 0) {
          second_result_bin[i] <- original_cumsum_new[i]/denominator_bin[i]
        } else{
          second_result_bin[i] <- 0
        }
      } 
      
      second_result_bin[second_result_bin == "NaN"] <- 0 
      second_result_bin[second_result_bin == "Inf"] <- 0 
      sd_second_bin[is.na(sd_second_bin)] <- 0
      
      first_non_zero  <- which(deg_bin_new > 0)[1]
      while (second_result_bin[first_non_zero] == 0) first_non_zero <- first_non_zero + 1
      factor_second <- second_result_bin[first_non_zero]
      second_result_bin <- second_result_bin / second_result_bin[first_non_zero]
      
      
      
      second_result_bin[second_result_bin == "NaN"] <- 0 
      second_result_bin[second_result_bin == "Inf"] <- 0 
      
      sd_result_second_bin <- original_cumsum_new * sd_second_bin/denominator_bin^2 / factor_second
      sd_result_second_bin[sd_result_second_bin == "NaN"] <- 0 
      sd_result_second_bin[sd_result_second_bin == "Inf"] <- 0 
      sd_log_result_second_bin <- sd_result_second_bin/second_result_bin
      sd_log_result_second_bin[sd_log_result_second_bin == "NaN"] <- 0 
      sd_log_result_second_bin[sd_log_result_second_bin == "Inf"] <- 0 
      upper_bin <- exp(log(second_result_bin) + 2 * sd_log_result_second_bin)
      lower_bin <- exp(log(second_result_bin) - 2 * sd_log_result_second_bin)
      
      
      
      regress_res_bin <- .regress_alpha(k = bin_object_new$center_bin, 
                                       A = second_result_bin,
                                       sd_log_A = sd_log_result_second_bin)
      alpha_second_bin <- regress_res_bin$alpha
      
      second_result_with_bin <- list(k = bin_object_new$center_bin,
                                     A = second_result_bin,
                                     A_upper = upper_bin,
                                     A_lower = lower_bin,
                                     alpha = alpha_second_bin,
                                     regress_res = regress_res_bin)
    }
  }
  if ("SI" == type) {
    # selective inference
    {
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
          var(bin_dist[i] * final_appear[,deg[i] + 1])/M
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
      
      
      selective_result_with_bin <- list(A = selective_result_bin,
                                        k = bin_object_new$center_bin,
                                        alpha = alpha_select_bin,
                                        alpha_no_sd = alpha_select_bin_no_sd,
                                        A_upper = upper,
                                        A_lower = lower,
                                        regress_result = regress_select_bin,
                                        regress_result_no_sd = regress_select_bin_no_sd)
      
    }
  }
  
  ############################################################
  ############################################################
  
  
  
  return(list(second_result_with_bin = second_result_with_bin,
              selective_result_with_bin = selective_result_with_bin,
              g_new = g_new,
              g = g,
              bin_object_old = bin_object,
              bin_object = bin_object_new, 
              appear_second = appear_second,
              final_deg = bin_dist,
              final_appear = final_appear,
              deg = deg,
              type = type,
              M = M,
              final_PA_type = final_PA_type,
              sd_selective_bin = sd_selective_bin,
              sd_result_second_bin = sd_result_second_bin))
  
}