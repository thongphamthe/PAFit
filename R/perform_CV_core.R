.performCV_core <- function(cv_data                                                  ,
                            r              = 10^c(-5 , -4 , -3 , -2, - 1, 0)       ,
                            s              = 10^c(-2 , -1 ,  0 ,  1 ,  2 ,  3 ,  4)  , 
                            stop_cond      = 10^-6                                   ,
                            print_out      = FALSE                                   ,
                            cv_deg_thresh  = c(1,10)                                 ,
                            normal_start_f = TRUE                                    ,
                            weight_f       = 0                                       ,
                           ...) { 
  
  FitMultinomial         <- function(true,dat){
    true[true == 0] <- 1
    return(sum(dat*log10(true)))
  }
 
  
  count <- 0
  
  #### Three passes to find optimal s #######
  #### for each of the cv_deg_thresh  #######
  rate_PAFit             <- s
  
  alpha_each           <- matrix(0,nrow = length(cv_deg_thresh), ncol = length(rate_PAFit))
  colnames(alpha_each) <- rate_PAFit
  
  estimated_fitness    <- list(length = length(cv_deg_thresh))
  #estimated_PA        <- NULL
  max_val              <- rep(-Inf,length(cv_deg_thresh))
  alpha_optimal        <- rep(-1,length(cv_deg_thresh))
  s_optimal            <- rep(0,length(cv_deg_thresh))
  

  
  for (j in 1:length(rate_PAFit)) {
    
    #print(paste0("Processing case ",count, " of ",total))
    result_PAFit  <- PAFit(cv_data$stats     , mode_f    = "Log_linear"    , s = rate_PAFit[j], 
                           auto_stop      =  TRUE , stop_cond  = stop_cond , normalized_f = FALSE ,
                           normal_start_f = normal_start_f,
                           weight_power   = weight_f,
                           ...);

    
    
    alpha_temp    <- result_PAFit$alpha
    #print(alpha_temp)
    if (TRUE == result_PAFit$diverge_zero) {
      alpha_each[,j] <- -Inf
    }
    else  { 
        for (dd in 1:length(cv_deg_thresh)) { 
            alpha_each[dd,j] <- 0  
            chosen_node      <- names(cv_data$stats$z_j[cv_data$stats$z_j >= cv_deg_thresh[dd]])
            for (k in 1:length(cv_data$m_each)) 
                if (cv_data$m_each[k] != 0) {
                    #PA              <- result_PAFit$A[cv_data$deg_each[k,] + 1] 
                    #PA[is.na(PA)]   <- result_PAFit$A[length(result_PAFit$A)]
                    PA              <-  cv_data$deg_each[k,chosen_node]^alpha_temp
                    #print(PA)
                    PA[PA == 0]     <- 1
                    fitness         <- rep(1,dim(cv_data$deg_each[,chosen_node, drop = FALSE])[2])
                    names(fitness)  <- colnames(cv_data$deg_each[,chosen_node, drop = FALSE])
                    fitness[chosen_node] <- result_PAFit$f[chosen_node] 
        
                    prob_PAFit      <- PA * fitness
                    prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
                    prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
                    alpha_each[dd,j]   <- alpha_each[dd, j] + 
                                           FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,chosen_node] * 
                                                          cv_data$m_each[k])) 
      }
        }
    }
    for (dd in 1:length(cv_deg_thresh)) 
    if (alpha_each[dd,j] > max_val[dd]) {
      estimated_fitness[[dd]]  <- result_PAFit$f[as.character(cv_data$stats$f_position)]
      s_optimal[dd]            <- rate_PAFit[j]
      max_val[dd]              <- alpha_each[dd,j]
      alpha_optimal[dd]        <- result_PAFit$alpha
    }
  }
  
  if (print_out == TRUE)
    print(paste0("After first iteration: s_optimal = ",s_optimal, " , alpha_optimal = ", alpha_optimal))  
  
  
  
  #s_index    <- which.max(alpha_each)[1]
  #s_optimal  <- rate_PAFit[s_index]
  for (dd in 1:length(cv_deg_thresh)) {
      chosen_node      <- names(cv_data$stats$z_j[cv_data$stats$z_j >= cv_deg_thresh[dd]]) 
      rate_PAFit <- c(0.2 , 0.4, 2.5, 5) * s_optimal[dd]
  
      #print(rate_PAFit)
  
      rate_PAFit <- rate_PAFit[rate_PAFit <= 10^4]
  
      rate_PAFit <- rate_PAFit[rate_PAFit >= 10^-2]

  
      alpha_each           <- rep(0,length(rate_PAFit))
      names(alpha_each)    <- rate_PAFit
  
      #estimated_PA        <- NULL
  
  #max_val    <- -Inf
  if (length(rate_PAFit) > 0)
    for (j in 1:length(rate_PAFit)) {
      #print(paste0("Processing case ",count, " of a maximum of ",total))
      result_PAFit <- PAFit(cv_data$stats, mode_f = "Log_linear", 
                            s = rate_PAFit[j], 
                            start_f      = estimated_fitness[[dd]], 
                            alpha_start  = alpha_optimal[dd], 
                            auto_stop    = TRUE, 
                            stop_cond    = stop_cond, 
                            normalized_f = FALSE, 
                            normal_start_f = normal_start_f,
                            weight_power   = weight_f,
                            ...);
      alpha_temp    <- result_PAFit$alpha  
      if (TRUE == result_PAFit$diverge_zero)
        alpha_each[j] <- -Inf
      else { 
          for (k in 1:length(cv_data$m_each)) 
              if (cv_data$m_each[k] != 0) {
              #PA              <- result_PAFit$A[cv_data$deg_each[k,] + 1] 
              #PA[is.na(PA)]   <- result_PAFit$A[length(result_PAFit$A)]
              PA              <-  cv_data$deg_each[k,chosen_node]^alpha_temp
              #print(PA)
              PA[PA == 0]     <- 1
              fitness         <- rep(1,dim(cv_data$deg_each[,chosen_node, drop = FALSE])[2])
              names(fitness)  <- colnames(cv_data$deg_each[,chosen_node, drop = FALSE])
              fitness[chosen_node] <- result_PAFit$f[chosen_node] 
              prob_PAFit      <- PA * fitness
              prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
              prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
              alpha_each[j] <- alpha_each[j] + 
                               FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,chosen_node] * 
                                                                           cv_data$m_each[k]))  
              }
        }
      if (alpha_each[j] > max_val[dd]) {
        estimated_fitness[[dd]] <- result_PAFit$f[as.character(cv_data$stats$f_position)]
        s_optimal[dd]           <- rate_PAFit[j]
        max_val[dd]             <- alpha_each[j]
        alpha_optimal[dd]       <- result_PAFit$alpha
      }
    }
  }
  #print(alpha_each)  
  
  if (print_out == TRUE)
    print(paste0("After second iteration: s_optimal = ",s_optimal, " , alpha_optimal = ", alpha_optimal))  
  
  
  for (dd in 1:length(cv_deg_thresh)) {
      rate_PAFit <- c(1/1.75 , 1/1.5 , 1/1.25 ,1, 1.25 , 1.5 , 1.75) * s_optimal[dd]
  
  
      rate_PAFit <- rate_PAFit[rate_PAFit <= 10^4]
  
      rate_PAFit <- rate_PAFit[rate_PAFit >= 10^-2]
  
  
       chosen_node      <- names(cv_data$stats$z_j[cv_data$stats$z_j >= cv_deg_thresh[dd]]) 
    
    #print(rate_PAFit)
    
    
    
    alpha_each           <- rep(0,length(rate_PAFit))
    names(alpha_each)    <- rate_PAFit
    
    #estimated_PA        <- NULL
    
    #max_val    <- -Inf
    if (length(rate_PAFit) > 0)
      for (j in 1:length(rate_PAFit)) {
        #print(paste0("Processing case ",count, " of a maximum of ",total))
        result_PAFit <- PAFit(cv_data$stats, mode_f = "Log_linear", 
                              s = rate_PAFit[j], 
                              start_f      = estimated_fitness[[dd]], 
                              alpha_start  = alpha_optimal[dd], 
                              auto_stop    = TRUE, 
                              stop_cond    = stop_cond, 
                              normal_start_f = normal_start_f,
                              weight_power   = weight_f ,
                              normalized_f = FALSE, ...);
        alpha_temp    <- result_PAFit$alpha  
        if (TRUE == result_PAFit$diverge_zero)
          alpha_each[j] <- -Inf
        else { 
          for (k in 1:length(cv_data$m_each)) 
            if (cv_data$m_each[k] != 0) {
              #PA              <- result_PAFit$A[cv_data$deg_each[k,] + 1] 
              #PA[is.na(PA)]   <- result_PAFit$A[length(result_PAFit$A)]
              PA              <-  cv_data$deg_each[k,chosen_node]^alpha_temp
              #print(PA)
              PA[PA == 0]     <- 1
              fitness         <- rep(1,dim(cv_data$deg_each[,chosen_node, drop = FALSE])[2])
              names(fitness)  <- colnames(cv_data$deg_each[,chosen_node, drop = FALSE])
              fitness[chosen_node] <- result_PAFit$f[chosen_node] 
              prob_PAFit      <- PA * fitness
              prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
              prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
              alpha_each[j] <- alpha_each[j] + 
                FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,chosen_node] * 
                                                                               cv_data$m_each[k]))  
            }
        }
        if (alpha_each[j] > max_val[dd]) {
          estimated_fitness[[dd]] <- result_PAFit$f[as.character(cv_data$stats$f_position)]
          s_optimal[dd]           <- rate_PAFit[j]
          max_val[dd]             <- alpha_each[j]
          alpha_optimal[dd]       <- result_PAFit$alpha
        }
      }
  }
  #print(alpha_each)  
  if (print_out == TRUE)
    print(paste0("Finally: s_optimal = ",s_optimal, " , alpha_optimal = ", alpha_optimal))  
  
  
  ok_index                <- 1:length(cv_deg_thresh)
  s_optimal_final         <- mean(s_optimal[ok_index])
  alpha_optimal_final     <- mean(alpha_optimal[ok_index])
  estimated_fitness_final <- estimated_fitness[[ok_index[1]]]
  if (length(ok_index) > 1)
    for (temp_jj in 2:length(ok_index))
      estimated_fitness_final <- estimated_fitness_final + estimated_fitness[[ok_index[temp_jj]]] 
  estimated_fitness_final <- estimated_fitness_final / length(ok_index)
  
  
  ### One pass to find optimal r #######
  chosen_node            <- names(cv_data$stats$z_j[cv_data$stats$z_j >= cv_deg_thresh[1]])
  
  ratio_vec_PAFit        <- sort(r,decreasing = TRUE)
  PA_each                <- rep(0,length(ratio_vec_PAFit))
  names(PA_each)         <- ratio_vec_PAFit
  max_val                <- -Inf 
  estimated_PA           <- NULL
  
  switch_flag            <- 0
  
  for (i in 1:length(ratio_vec_PAFit)) {
    count <- count + 1
    #print(paste0("Processing case ",count, " of a maximum of ",total))
    if (!is.null(estimated_PA))
      result_PAFit <- PAFit(cv_data$stats, 
                            s = s_optimal_final, 
                            r = ratio_vec_PAFit[i], 
                            only_PA   = TRUE, 
                            true_f    = estimated_fitness_final, 
                            auto_stop =  TRUE, 
                            start_A   = estimated_PA,
                            stop_cond = stop_cond, 
                            normal_start_f = normal_start_f,
                            weight_power   = weight_f ,
                            normalized_f = FALSE, ...)
    else result_PAFit <- PAFit(cv_data$stats, 
                          s = s_optimal_final, 
                          r = ratio_vec_PAFit[i], 
                          only_PA = TRUE, 
                          true_f = estimated_fitness_final, 
                          auto_stop =  TRUE, 
                          alpha_start = alpha_optimal_final,
                          stop_cond = stop_cond, 
                          normal_start_f = normal_start_f,
                          normalized_f = FALSE, ...)
    
    
    alpha_temp    <- result_PAFit$alpha  
    
    if (TRUE == result_PAFit$diverge_zero)
       PA_each[i] <- -Inf
    else for (k in 1:length(cv_data$m_each))
      if (cv_data$m_each[k] != 0) { 
        #PA              <- cv_data$deg_each[k,chosen_node]^alpha_temp
        PA              <- result_PAFit$A[cv_data$deg_each[k,chosen_node]]
        PA[is.na(PA)]   <- cv_data$deg_each[k,chosen_node][is.na(PA)]^alpha_temp
        #print(PA)
        PA[PA == 0]     <- mean(result_PAFit$A)
        fitness         <- rep(1,dim(cv_data$deg_each[,chosen_node, drop = FALSE])[2])
        names(fitness)  <- colnames(cv_data$deg_each[,chosen_node, drop = FALSE])
        fitness[chosen_node] <- result_PAFit$f[chosen_node] 
        
        prob_PAFit      <- PA * fitness
        prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
        prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
        PA_each[i]      <- PA_each[i] + 
          FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,chosen_node] * 
                                                                         cv_data$m_each[k])) 
      }
    if (i == 1) {  
        lambda_optimal    <- result_PAFit$lambda  
        r_optimal         <- ratio_vec_PAFit[i]
        max_val           <- PA_each[i]
        estimated_PA      <- result_PAFit$theta
    } else if (PA_each[i] < PA_each[i - 1]) {
        break  
    } else {
        lambda_optimal    <- result_PAFit$lambda  
        r_optimal         <- ratio_vec_PAFit[i]
        max_val           <- PA_each[i]
        estimated_PA      <- result_PAFit$theta  
    }
  }
  
  
  #r_index    <- which.max(PA_each)[1]
  #r_optimal  <- ratio_vec_PAFit[r_index]
  
  # ratio_vec_PAFit <- c(0.2 , 0.4 , 0.6 , 0.8 , 2 , 4 , 6 , 8) * r_optimal 
  # 
  # ratio_vec_PAFit <- ratio_vec_PAFit[ratio_vec_PAFit <= 0.1]
  # 
  # #ratio_vec_PAFit <- ratio_vec_PAFit[ratio_vec_PAFit >= 10^-5]
  # 
  # PA_each         <- rep(0,length(ratio_vec_PAFit))
  # names(PA_each)  <- ratio_vec_PAFit
  # 
  # if (length(ratio_vec_PAFit) > 0) {
  #   
  #   for (i in 1:length(ratio_vec_PAFit)) {
  #     count <- count + 1
  #     #print(paste0("Processing case ",count, " of a maximum of ",total))
  #     result_PAFit <- PAFit(cv_data$stats,
  #                           s           = s_optimal_final, 
  #                           r           = ratio_vec_PAFit[i], 
  #                           auto_stop   = TRUE, 
  #                           start_A     = estimated_PA,
  #                           true_f      = estimated_fitness_final,  
  #                           only_PA     = TRUE, 
  #                           #alpha_start = result_temp$alpha,
  #                           stop_cond   = stop_cond       ,
  #                           ...)
  #     alpha_temp    <- result_PAFit$alpha    
  #     if (TRUE == result_PAFit$diverge_zero)
  #       PA_each[j] <- -Inf
  #     else for (k in 1:length(cv_data$m_each))
  #       if (cv_data$m_each[k] != 0) { 
  #           PA              <- result_PAFit$A[cv_data$deg_each[k,chosen_node]]
  #           PA[is.na(PA)]   <- cv_data$deg_each[k,chosen_node][is.na(PA)]^alpha_temp
  #         fitness         <- rep(1,dim(cv_data$deg_each[,chosen_node])[2])
  #         names(fitness)  <- colnames(cv_data$deg_each[,chosen_node])
  #         fitness[chosen_node] <- result_PAFit$f[chosen_node] 
  #         
  #         prob_PAFit      <- PA * fitness
  #         prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
  #         prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
  #         alpha_each[j] <- alpha_each[j] + 
  #           FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,chosen_node] * 
  #                                                                          cv_data$m_each[k])) 
  #       }
  #     if (PA_each[i] > max_val) {
  #       r_optimal         <- ratio_vec_PAFit[i]
  #       lambda_optimal    <- result_PAFit$lambda  
  #       max_val           <- PA_each[i]
  #       estimated_PA      <- result_PAFit$theta
  #     }
  #   }
  # }
  result    <- list(r_optimal         = r_optimal, 
                    lambda_optimal    = lambda_optimal, 
                    s_optimal         = s_optimal_final,
                    alpha_optimal     = alpha_optimal_final,
                    estimated_fitness = estimated_fitness_final,
                    estimated_PA      = estimated_PA,
                    cv_deg_thresh     = cv_deg_thresh)
  class(result) <- "CV_Result"
  return(result)
}


