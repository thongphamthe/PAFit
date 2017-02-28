.performCV_new <- function(cv_data                                          ,
                       r              = 10^c(-5, -4 , -3 , -2 , -1 , 0 , 1) ,
                       s              = 10^c(-1 , 0 ,  1 ,  2 ,  3 , 4)     , 
                       stop_cond      = 10^-6                               ,
                       print_out      = FALSE                               ,                   
                       ...) { 
  
  FitMultinomial         <- function(true,dat){
    true[true == 0] <- 1
    return(sum(dat*log(true)))
  }
  
  
  count <- 0
  
  ### Three passes to find optimal s #######
  rate_PAFit             <- s
  alpha_each             <- rep(0,length(rate_PAFit))
  names(alpha_each)      <- rate_PAFit
  estimated_fitness      <- vector()
  estimated_PA           <- NULL
  max_val                <- -Inf 
  alpha_optimal          <- -1
  
  true_ll                <- alpha_optimal
  norm_ll                <- alpha_optimal  
  
  sum_using_true         <- alpha_optimal
  sum_using_norm         <- alpha_optimal
  
  
  for (j in 1:length(rate_PAFit)) {
    
    #print(paste0("Processing case ",count, " of ",total))
    result_PAFit  <- PAFit_new(cv_data$stats     , mode_f    = "Log_linear" , s = rate_PAFit[j], 
                           auto_stop =  TRUE , stop_cond  = stop_cond     , normalized_f = FALSE , ...);
    alpha_each[j] <- 0
    
    true_ll[j]    <-  result_PAFit$unnorm_ll 
    
    norm_ll[j]    <-  result_PAFit$norm_ll 
    
    alpha_temp    <- result_PAFit$alpha
    #print(alpha_temp)
    if (TRUE == result_PAFit$diverge_zero)
      alpha_each[j] <- -Inf
    else for (k in 1:length(cv_data$m_each)) 
      if (cv_data$m_each[k] != 0) {
        #PA              <- result_PAFit$A[cv_data$deg_each[k,] + 1] 
        #PA[is.na(PA)]   <- result_PAFit$A[length(result_PAFit$A)]
        PA              <-  cv_data$deg_each[k,]^alpha_temp
        #print(PA)
        PA[PA == 0]     <- 1
        fitness         <- rep(1,dim(cv_data$deg_each)[2])
        names(fitness)  <- colnames(cv_data$deg_each)
        fitness[as.character(cv_data$stats$f_position)] <- result_PAFit$f[as.character(cv_data$stats$f_position)] 
 
        prob_PAFit      <- PA * fitness
        prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
        prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
        alpha_each[j] <- alpha_each[j] + 
          FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,] * 
                                                                         cv_data$m_each[k])) 
      }
    #sum_using_true[j] <- true_ll[j]  + alpha_each[j]
    
    #sum_using_norm[j] <- norm_ll[j] + alpha_each[j]/sum(cv_data$m_each)
    
    #alpha_ea
    if (alpha_each[j] > max_val) {
      estimated_fitness <- result_PAFit$f[as.character(cv_data$stats$f_position)]
      s_optimal         <- rate_PAFit[j]
      max_val           <- alpha_each[j]
      alpha_optimal     <- result_PAFit$alpha
    }
  }
  #print(alpha_each)  
  
  #print("True ll:")
  #print(true_ll)
  
  #print("Norm ll:")
  #print(norm_ll)
  
  #print("Sum Using True ll:")
  #print(sum_using_true)
  
  
  #print("Sum using norm ll:")
  #print(sum_using_norm)
  
  #print(which.max(norm_ll))
  
  
  if (print_out == TRUE)
    print(paste0("After first iteration: s_optimal = ",s_optimal, " , alpha_optimal = ", alpha_optimal))  
  
  #s_index    <- which.max(alpha_each)[1]
  #s_optimal  <- rate_PAFit[s_index]
  rate_PAFit <- c(0.2 , 0.4 , 2.5, 5) * s_optimal
  
  #print(rate_PAFit)
  
  rate_PAFit <- rate_PAFit[rate_PAFit <= 10^4]
  
  rate_PAFit <- rate_PAFit[rate_PAFit >= 10^-1]
  
  alpha_each             <- rep(0,length(rate_PAFit))
  names(alpha_each)      <- rate_PAFit
  #max_val    <- -Inf
  if (length(rate_PAFit) > 0)
    for (j in 1:length(rate_PAFit)) {
      #print(paste0("Processing case ",count, " of a maximum of ",total))
      result_PAFit <- PAFit_new(cv_data$stats, mode_f = "Log_linear", s = rate_PAFit[j], 
                            start_f      = estimated_fitness, 
                            alpha_start  = alpha_optimal, 
                            auto_stop    = TRUE, 
                            stop_cond    = stop_cond, 
                            normalized_f = FALSE,...);
      alpha_each[j] <- 0;
      alpha_temp    <- result_PAFit$alpha  
      if (TRUE == result_PAFit$diverge_zero)
        alpha_each[j] <- -Inf
      else for (k in 1:length(cv_data$m_each)) 
        if (cv_data$m_each[k] != 0) {
          #PA              <- result_PAFit$A[cv_data$deg_each[k,] + 1] 
          #PA[is.na(PA)]   <- result_PAFit$A[length(result_PAFit$A)]
          PA              <-  cv_data$deg_each[k,]^alpha_temp
          #print(PA)
          PA[PA == 0]     <- 1
          fitness         <- rep(1,dim(cv_data$deg_each)[2])
          names(fitness)  <- colnames(cv_data$deg_each)
          fitness[as.character(cv_data$stats$f_position)] <- result_PAFit$f[as.character(cv_data$stats$f_position)] 
          
          prob_PAFit      <- PA * fitness
          
          
          prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
          prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
          alpha_each[j]   <- alpha_each[j] + 
            FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,] * cv_data$m_each[k])) 
        }
      if (alpha_each[j] > max_val) {
        estimated_fitness <- result_PAFit$f[as.character(cv_data$stats$f_position)]
        s_optimal         <- rate_PAFit[j]
        max_val           <- alpha_each[j]
        alpha_optimal     <- result_PAFit$alpha
      }
    }
  #print(alpha_each)  
  
  if (print_out == TRUE)
    print(paste0("After second iteration: s_optimal = ",s_optimal, " , alpha_optimal = ", alpha_optimal))  
  
  
  rate_PAFit <- c(1/1.75 , 1/1.5 , 1/1.25 , 1.25 , 1.5 , 1.75) * s_optimal
  
  
  rate_PAFit <- rate_PAFit[rate_PAFit <= 10^4]
  
  rate_PAFit <- rate_PAFit[rate_PAFit >= 10^-1]
  
  alpha_each             <- rep(0,length(rate_PAFit))
  names(alpha_each)      <- rate_PAFit
  
  if (length(rate_PAFit) > 0)
    for (j in 1:length(rate_PAFit)) {
      #print(paste0("Processing case ",count, " of a maximum of ",total))
      result_PAFit <- PAFit_new(cv_data$stats, mode_f = "Log_linear", s = rate_PAFit[j], 
                                start_f = estimated_fitness, 
                                alpha_start = alpha_optimal, 
                                auto_stop =  TRUE, stop_cond = stop_cond, normalized_f = FALSE,...);
      alpha_each[j] <- 0;
      alpha_temp    <- result_PAFit$alpha
      if (TRUE == result_PAFit$diverge_zero)
        alpha_each[j] <- -Inf
      else for (k in 1:length(cv_data$m_each)) 
        if (cv_data$m_each[k] != 0) {
          #PA              <- result_PAFit$A[cv_data$deg_each[k,] + 1] 
          #PA[is.na(PA)]   <- result_PAFit$A[length(result_PAFit$A)]
          PA              <-  cv_data$deg_each[k,]^alpha_temp
          #print(PA)
          PA[PA == 0]     <- 1
          fitness         <- rep(1,dim(cv_data$deg_each)[2])
          names(fitness)  <- colnames(cv_data$deg_each)
          fitness[as.character(cv_data$stats$f_position)] <- result_PAFit$f[as.character(cv_data$stats$f_position)] 
          
          prob_PAFit      <- PA * fitness
          
          
          prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
          prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
          alpha_each[j] <- alpha_each[j] + 
            FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,] * cv_data$m_each[k])) 
        }
      
      if (alpha_each[j] > max_val) {
        estimated_fitness <- result_PAFit$f[as.character(cv_data$stats$f_position)]
        s_optimal         <- rate_PAFit[j]
        max_val           <- alpha_each[j]
        alpha_optimal     <- result_PAFit$alpha
      }
    }
  #print(alpha_each)  
  if (print_out == TRUE)
    print(paste0("Finally: s_optimal = ",s_optimal, " , alpha_optimal = ", alpha_optimal))  
  
  #print(paste0("s_optimal: ", s_optimal))
  #print(alpha_optimal)
  #s_index    <- which.max(alpha_each)[1]
  #s_optimal  <- rate_PAFit[s_index]
  
  
  ### Two passes to find optimal r #######
  ratio_vec_PAFit        <- r
  PA_each                <- rep(0,length(ratio_vec_PAFit))
  names(PA_each)         <- ratio_vec_PAFit
  max_val                <- -Inf 
  for (i in 1:length(ratio_vec_PAFit)) {
    count <- count + 1
    #print(paste0("Processing case ",count, " of a maximum of ",total))
    result_PAFit <- PAFit_new(cv_data$stats, 
                          s = s_optimal, 
                          r = ratio_vec_PAFit[i], 
                          only_PA = TRUE, 
                          true_f = estimated_fitness, 
                          auto_stop =  TRUE, 
                          alpha_start = alpha_optimal,
                          stop_cond = stop_cond, 
                          normalized_f = FALSE,...)
    
    
    alpha_temp    <- result_PAFit$alpha  
    
    if (TRUE == result_PAFit$diverge_zero)
      PA_each[j] <- -Inf
    else for (k in 1:length(cv_data$m_each))
      if (cv_data$m_each[k] != 0) { 
        PA              <-  cv_data$deg_each[k,]^alpha_temp
        #print(PA)
        PA[PA == 0]     <- 1
        
        fitness         <- rep(1,dim(cv_data$deg_each)[2])
        names(fitness)  <- colnames(cv_data$deg_each)
        
        fitness[as.character(cv_data$stats$f_position)] <- result_PAFit$f[as.character(cv_data$stats$f_position)] 
        
        prob_PAFit      <- PA * fitness
        
        
        
        prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE); 
        prob_PAFit[sapply(prob_PAFit,is.na)] <- 0; 
        PA_each[i]      <- PA_each[i] + 
          FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,] * cv_data$m_each[k])) ;
      }
    if (PA_each[i] > max_val) {
      r_optimal         <- ratio_vec_PAFit[i]
      max_val           <- PA_each[i]
      estimated_PA      <- result_PAFit$theta
    }    
  }
  
  
  #r_index    <- which.max(PA_each)[1]
  #r_optimal  <- ratio_vec_PAFit[r_index]
  
  ratio_vec_PAFit <- c(0.2 , 0.4 , 0.6 , 0.8 , 2 , 4 , 6 , 8) * r_optimal 
  
  ratio_vec_PAFit <- ratio_vec_PAFit[ratio_vec_PAFit <= 10]
  
  PA_each         <- rep(0,length(ratio_vec_PAFit))
  names(PA_each)  <- ratio_vec_PAFit
  
  if (length(ratio_vec_PAFit) > 0) {
    
    for (i in 1:length(ratio_vec_PAFit)) {
      count <- count + 1
      #print(paste0("Processing case ",count, " of a maximum of ",total))
      result_PAFit <- PAFit_new(cv_data$stats,
                                s           = s_optimal, 
                                r           = ratio_vec_PAFit[i], 
                                auto_stop   = TRUE, 
                                start_A     = estimated_PA,
                                true_f      = estimated_fitness ,  
                                only_PA     = TRUE, 
                               #alpha_start = result_temp$alpha,
                                stop_cond   = stop_cond        ,
                                ...)
      alpha_temp    <- result_PAFit$alpha    
      if (TRUE == result_PAFit$diverge_zero)
        PA_each[j] <- -Inf
      else for (k in 1:length(cv_data$m_each))
        if (cv_data$m_each[k] != 0) { 
          PA              <-  cv_data$deg_each[k,]^alpha_temp
          #print(PA)
          PA[PA == 0]     <- 1
          
          fitness         <- rep(1,dim(cv_data$deg_each)[2])
          names(fitness)  <- colnames(cv_data$deg_each)
          
          fitness[as.character(cv_data$stats$f_position)] <- result_PAFit$f[as.character(cv_data$stats$f_position)] 
          
          prob_PAFit      <- PA * fitness
          
          
          prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE); 
          prob_PAFit[sapply(prob_PAFit,is.na)] <- 0; 
          PA_each[i]    <- PA_each[i] + 
            FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,] * cv_data$m_each[k])) ;
        }
      if (PA_each[i] > max_val) {
        r_optimal         <- ratio_vec_PAFit[i]
        max_val           <- PA_each[i]
        estimated_PA      <- result_PAFit$theta
      }
    }
  }
  result    <- list(r_optimal         = r_optimal, 
                    alpha_optimal     = alpha_optimal,
                    s_optimal         = s_optimal,
                    estimated_fitness = estimated_fitness,
                    estimated_PA      = estimated_PA)
  class(result) <- "CV_Result"
  #print(paste0("Optimal r parameter is: ",r_optimal));
  #print(paste0("Optimal s parameter is: ",s_optimal));
  
  return(result)
}


