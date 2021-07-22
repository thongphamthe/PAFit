
.estimate_new_correct_mcPAMLE <- function(net, M = 5,
                         net_type = "directed",
                         true_time  = NULL,
                         final_PA_type = "mean",
                         G = 100) {
  
  #result_use = "Iterative"
  #which_interpolate = "GeoMean" #, "LogLinear"
  
  if (M <= 2)
    stop(" M should be at least 5.")
  {
    from_node      <- net[,1]
    to_node        <- net[net[,2] != -1,2]
    unique_node    <- unique(c(from_node,to_node))
    deg_vec        <- rep(0,length(unique_node))
    names(deg_vec) <- as.character(as.integer(unique_node))
    if (net_type == "directed") {
      hist           <- table(to_node)
    } else if (net_type == "undirected") {
      #get self-loop
      edge     <- net[net[,2] != -1,1:2,drop = FALSE]
      self     <- edge[edge[,1] == edge[,2],1]
      non_self <-  edge[edge[,1] != edge[,2],1:2,drop = FALSE]
      # only count self-loop once
      hist <- table(c(self,non_self[,1],non_self[,2]))
    } else {
      stop("Wrong net_type")
    }
    deg_vec[names(hist)] <- hist
    Total_T        <- dim(net)[1] - 2 # remove the first time-step
    #print(Total_T)
    total_edge     <- sum(hist) 
    p_estimate     <- ((length(unique_node) - 2)) / ((length(unique_node) - 2) + (total_edge - 1))
    
    #print(p_estimate)
    
    
    dist         <- table(deg_vec)
    deg_max      <- max(as.integer(names(dist)))
    deg          <- as.integer(names(dist[-length(dist)]))
    names(deg)   <- as.integer(deg)
    
    deg_penmax   <- max(deg)
    
    
    
    }
  
  
    
    # first do it with no binning
    g              <- deg_penmax + 1  # from 0 to deg_penmax
    bin_object     <- .binning(deg_penmax,G = g)
    g              <- bin_object$G_small
    
    degree_exist           <- rep(0,g)
    names(degree_exist)    <- as.integer(0:(g-1))
    temp_dist              <- dist[-length(dist)]
    
    degree_exist[temp_dist[names(degree_exist)] > 0] <- 1
    
    
    
    bin_dist       <- rep(0,g)
    for (i in 1:g){
      index <- which(deg >= bin_object$start[i] & deg <= bin_object$end[i])
      bin_dist[i] <- sum(dist[index])
    }
    names(bin_dist) <- 0:(g-1)
    deg_bin         <- bin_object$center_bin
    names(deg_bin)  <- 0:(g-1)
    
    #first estimate
    first_result        <- rep(1,g)
    names(first_result) <- 0:(g - 1) 
    
    E       <- sum(deg_vec)
    
    original_cumsum <- rep(0,length(first_result))
    names(original_cumsum) <- 0:(g - 1) 
    
    for (i in 1:(length(dist) - 1)) {
      temp <- sum(dist[(i + 1):length(dist)]) # the number of new edges connecting to that degree i 
      bin_index <- which(bin_object$start <= deg[i] & deg[i] <= bin_object$end) # find the bin of that degree i
      if (length(bin_index) != 1) print("Wrong at bin index")
      original_cumsum[bin_index] <- original_cumsum[bin_index] + temp;    # we have to add temp to this bin
    }
    

    #W_k     <- E - deg_bin
    W_k     <- rep(1,length(deg_bin))

    for ( gg in 1:g) {
      first_result[gg]    <- original_cumsum[gg]/bin_dist[gg]/W_k[gg]
    }
    first_appear <- W_k
    first_appear <- first_appear/first_appear[2]
    # calculate the silver result using true appearance time when there is no binning
    silver_result <- NULL
    if (g == deg_penmax + 1) {
      if (!is.null(true_time)) {
        silver_result <- original_cumsum/bin_dist/true_time[as.character(deg_bin)]
        silver_result[silver_result == "NaN"] <- 0 
        silver_result[silver_result == "Inf"] <- 0 
        
        
        first_non_zero <- 1
        while (silver_result[first_non_zero] == 0) first_non_zero <- first_non_zero + 1
        silver_result <- silver_result / silver_result[first_non_zero]
        names(silver_result) <- deg_bin
        silver_result <- list(k = deg_bin, A = silver_result)
      }
    }
    
    first_result[first_result == "NaN"] <- 0 
    first_result[first_result == "Inf"] <- 0 
    
    #first_non_zero <- 1  
    first_non_zero  <- which(deg_bin == 1)
    while (first_result[first_non_zero] == 0) first_non_zero <- first_non_zero + 1
    first_result <- first_result / first_result[first_non_zero]
    
    
    first_result[first_result == "NaN"] <- 0 
    first_result[first_result == "Inf"] <- 0 
    
  M2 <- M
  
  
  is_directed <- 1
  if (net_type == "directed") {
    is_directed <- 1
  } else if (net_type == "undirected") {
    is_directed <- 0
  } else {
    stop("Wrong")
  }
  #print(paste0("Is_directed: ",is_directed))
  
  
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
          # if(first_result[before] < first_result[after]) {
          #     custom_PA[i] <-  first_result[after] * (1 - (1 - (first_result[before] / first_result[after])) * 
          #                                                (bin_object$center_bin[bin_index]/(bin_object$center_bin[before] + 1)) / (bin_object$center_bin[after] / (bin_object$center_bin[before] + 1)));
          # } else {
          # #print("in here")
          #     custom_PA[i] <-  first_result[before] * (1 - (1 - first_result[after] / first_result[before]) * 
          #         (bin_object$center_bin[bin_index]/(bin_object$center_bin[before] + 1)) / (bin_object$center_bin[after] / (bin_object$center_bin[before] + 1))) ;  
          # }
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
  
  first_non_zero  <- which(deg_bin == 1)
  while (custom_PA[first_non_zero] == 0) first_non_zero <- first_non_zero + 1;
  custom_PA <- custom_PA / custom_PA[first_non_zero]

  #prepare bin object
  {
    bin_object_new     <- .binning(deg_penmax,G = G)
    g_new              <- bin_object_new$G_small
    bin_dist_new       <- rep(0,g_new)
    names(bin_dist_new) <- 0:(g_new-1)
    deg_bin_new         <- bin_object_new$center_bin
    names(deg_bin_new)  <- 0:(g_new-1)
    original_cumsum_new <- rep(0,g_new)
    names(original_cumsum_new) <- 0:(g_new - 1) 
    for (i in 1:g_new){
      index <- which(deg >= bin_object_new$start[i] & deg <= bin_object_new$end[i])
      bin_dist_new[i] <- sum(dist[index])
    }
    for (i in 1:(length(dist) - 1)) {
      temp <- sum(dist[(i + 1):length(dist)]) # the number of new edges connecting to that degree i 
      bin_index_new <- which(bin_object_new$start <= deg[i] & deg[i] <= bin_object_new$end) # find the bin of that degree i
      if (length(bin_index_new) != 1) print("Wrong at bin index")
      original_cumsum_new[bin_index_new] <- original_cumsum_new[bin_index_new] + temp;    # we have to add temp to this bin
    }
  }
  
  appear_second <- matrix(0,nrow = M2, ncol = g)
  colnames(appear_second) <- 0:(g-1)
  final_appear <- matrix(0,nrow = M2, ncol = g)
  colnames(final_appear) <- 0:(g-1)
  
  #final_PA_val <- mean(custom_PA[custom_PA > 0])  
  #final_PA_val <- max(custom_PA[custom_PA > 0])  
  if ("mean" == final_PA_type) {
    final_PA_val <- mean(custom_PA[custom_PA > 0])  
  } else if ("max" == final_PA_type) {
    final_PA_val <- max(custom_PA[custom_PA > 0])  
  }
  mode_PA    = 3;
  bin_vector = bin_object$bin
  alpha = 1;
  beta =  2;
  #print(g)
  .generate_net_C_with_count_multi_corrected(appear_second,
                                   final_appear,
                                   Total_T, mode_PA, alpha,beta, custom_PA, p_estimate,
                                   bin_vector , g, 
                                   deg_penmax,final_PA_val,is_directed,M,degree_exist)
  
  ############################################################
  ################### binning result #########################
  ############################################################
  

  
  #### baseline ######
  {
    W_k_new                 <- rep(0,g_new)
    names(W_k_new)          <- 0:(g_new-1)
    first_result_bin        <- rep(0,g_new)
    names(first_result_bin) <- 0:(g_new-1)
    
    for (i in 1:(length(dist) - 1)) {
      bin_index_new <- which(bin_object_new$start <= deg[i] & deg[i] <= bin_object_new$end) # find the bin of that degree i
      if (length(bin_index_new) != 1) print("Wrong at bin index")
      W_k_new[bin_index_new] <- W_k_new[bin_index_new] +bin_dist[i] * W_k[deg[i] + 1];    # we have to add temp to this bin
    }
    
    for (i in 1:g_new){
      first_result_bin[i] <- original_cumsum_new[i]/W_k_new[i]
    }
    
    first_result_bin[first_result_bin == "NaN"] <- 0 
    first_result_bin[first_result_bin == "Inf"] <- 0 
    
    #first_non_zero <- 1  
    first_non_zero  <- which(deg_bin_new > 0)[1]
    while (first_result_bin[first_non_zero] == 0) first_non_zero <- first_non_zero + 1
    first_result_bin <- first_result_bin / first_result_bin[first_non_zero]
    
    first_result_bin[first_result_bin == "NaN"] <- 0 
    first_result_bin[first_result_bin == "Inf"] <- 0 
    
    regress_first_bin      <- .regress_alpha(k = bin_object_new$center_bin, 
                                            A = first_result_bin)
    alpha_first_bin        <- regress_first_bin$alpha
    
    first_estimate_bin <- list(A = first_result_bin,
                               k = bin_object_new$center_bin,
                               alpha = alpha_first_bin,
                               A_upper = first_result_bin,
                               A_lower = first_result_bin,
                               regress_result = regress_first_bin)
  }
  
  
  
  ### proposed method ####
  {
    second_result_bin        <- rep(0,g_new)
    names(second_result_bin) <- 0:(g_new-1)
    
    
    denominator_bin          <- rep(0,g_new)
    names(denominator_bin)   <- 0:(g_new - 1) 
    
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
  
  
  # selective inference
  {
    select_appear_bin        <- rep(0,g_new)
    names(select_appear_bin) <-  0:(g_new - 1) 
    denominator_select_bin   <- rep(0,g_new)
    sd_selective_bin            <- rep(0,g_new)
    names(sd_selective_bin)     <- 0:(g_new - 1) 
    var_vec_select              <- rep(0,g_new)
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
  
  ############################################################
  ############################################################
  
  
  return(list(first_estimate_bin = first_estimate_bin,
              second_result_with_bin = second_result_with_bin,
              selective_result_with_bin = selective_result_with_bin,
              g_new = g_new,
              g = g,
              W_k = W_k,
              bin_object_old = bin_object,
              bin_object = bin_object_new, 
              appear_second = appear_second,
              final_deg = bin_dist,
              deg_vec = deg_vec,
              dist = dist,
              bin_dist = bin_dist,
              bin_dist_new = bin_dist_new,
              silver_result = silver_result,
              original_cumsum = original_cumsum,
              original_cumsum_bin = original_cumsum_new,
              
              first_appear = first_appear,
              deg_penmax = deg_penmax,
              final_appear = final_appear,
              custom_PA = custom_PA,
              bin_object_new = bin_object_new,
              deg = deg,
              M = M,
              Total_T = Total_T,
              deg_bin = deg_bin,
              deg_bin_new = deg_bin_new,
              net_type = net_type,
              sd_selective_bin = sd_selective_bin,
              sd_result_second_bin = sd_result_second_bin,
              p_estimate = p_estimate,
              deg_max = deg_max,
              degree_exist = degree_exist,
              use_var = use_var))
  
}