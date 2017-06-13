.OnlyA_CV <- function(cv_data                                   ,
                     r         = c(0,10^c(-6, -5, -4, -3,-2,-1)),
                     stop_cond = 10^-6                          ,
                     print_out = FALSE                          ,     
                     rough     = TRUE                           ,  
                    ...) { 
  
  FitMultinomial         <- function(true,dat){
    true[true == 0] <- 1
    return(sum(dat*log(true)))
  }
  count <- 0
  
  ### Two passes to find optimal r #######
  ratio_vec_PAFit        <- sort(r,decreasing = TRUE)
  PA_each                <- rep(0,length(ratio_vec_PAFit))
  names(PA_each)         <- ratio_vec_PAFit
  max_val                <- -Inf 
  # A quick run to estimate the attachment exponent
  #print("Reached here")
  result_temp  <-  PAFit(cv_data$stats,
                         mode_f       = "Log_linear",  
                         only_PA      = TRUE, 
                         alpha_start  = 0.5,
                         stop_cond    = stop_cond,
                         ...)
  #print("Not yet reached here")
  
  for (i in 1:length(ratio_vec_PAFit)) {
    count <- count + 1
    # A more refined run with warm-start
    result_PAFit <- PAFit(cv_data$stats, 
                          r            = ratio_vec_PAFit[i], 
                          only_PA      = TRUE, 
                          alpha_start  = result_temp$alpha,
                          stop_cond    = stop_cond,
                          ...)
    
    #print(length(result_PAFit$theta))
    #print(dim(cv_data$n_tk_each))
    l_theta         <- length(result_PAFit$theta)
    estimated_alpha <- result_PAFit$alpha
    l_n_tk          <- dim(cv_data$n_tk_each)[2]
    PA_term         <- result_PAFit$theta
    if (l_n_tk > l_theta)
        PA_term <- c(PA_term,cv_data$center_k[(l_theta + 1):l_n_tk]^estimated_alpha)  
    if (TRUE == result_PAFit$diverge_zero)
      PA_each[i] <- -Inf
    else for (k in 1:length(cv_data$m_each))
      if (cv_data$m_each[k] != 0) { 
        prob_PAFit      <- PA_term  * cv_data$n_tk_each[k,] 
        prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE); 
        prob_PAFit[sapply(prob_PAFit,is.na)] <- 0; 
        PA_each[i]      <- PA_each[i] + 
          FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,] * cv_data$m_each[k])) ;
      }
    if (i == 1) {
        r_optimal         <- ratio_vec_PAFit[i]
        max_val           <- PA_each[i]
        estimated_PA      <- result_PAFit$theta
    } else if (PA_each[i] < PA_each[i - 1])    {
        break  
    } else {
        r_optimal         <- ratio_vec_PAFit[i]
        max_val           <- PA_each[i]
        estimated_PA      <- result_PAFit$theta  
    }
  }
  #r_index    <- which.max(PA_each)[1]
  #r_optimal  <- ratio_vec_PAFit[r_index]
  if (rough == FALSE) {
      ratio_vec_PAFit <- sort(c(0.2 , 0.5, 2 , 5) * r_optimal,decreasing = TRUE)
  
      ratio_vec_PAFit <- ratio_vec_PAFit[ratio_vec_PAFit <= 10^-1]
  
      PA_each         <- rep(0,length(ratio_vec_PAFit))
      names(PA_each)  <- ratio_vec_PAFit
  
      if (length(ratio_vec_PAFit) > 0) {
    
      for (i in 1:length(ratio_vec_PAFit)) {
          count <- count + 1
          #print(paste0("Processing case ",count, " of a maximum of ",total))
          result_PAFit <- PAFit(cv_data$stats, 
                                r            = ratio_vec_PAFit[i], 
                                only_PA      = TRUE, 
                                start_A      = estimated_PA,
                                stop_cond    = stop_cond,
                                ...)
         l_theta         <- length(result_PAFit$theta)
         estimated_alpha <- result_PAFit$alpha
         l_n_tk          <- dim(cv_data$n_tk_each)[2]
         PA_term         <- result_PAFit$theta
         if (l_n_tk > l_theta)
            PA_term <- c(PA_term,cv_data$center_k[(l_theta + 1):l_n_tk]^estimated_alpha)  
        if (TRUE == result_PAFit$diverge_zero)
            PA_each[i] <- -Inf
        else for (k in 1:length(cv_data$m_each))
            if (cv_data$m_each[k] != 0) { 
                prob_PAFit      <- PA_term  * cv_data$n_tk_each[k,] 
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
  }
  }
  result    <- list(r_optimal         = r_optimal, 
                    alpha_optimal     = result_temp$alpha,
                    estimated_PA      = estimated_PA)
  class(result) <- "CV_Result"
  #print(paste0("Optimal r parameter is: ",r_optimal));
  #print(paste0("Optimal s parameter is: ",s_optimal));
  
  return(result)
}


