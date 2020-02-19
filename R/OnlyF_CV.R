.OnlyF_CV <- function(cv_data                                             ,
                     s              = 10^c(-1 , 0 ,  1 ,  2 ,  3 , 4)     , 
                     stop_cond      = 10^-5                               ,
                     print_out      = FALSE                               ,
                     alpha                                                ,
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
  max_val                <- -Inf 
  
  
  
  
  for (j in 1:length(rate_PAFit)) {
    
    #print(paste0("Processing case ",count, " of ",total))
    result_PAFit  <- PAFit(cv_data$stats     , mode_f    = "Log_linear" , 
                           s = rate_PAFit[j], 
                           only_f = TRUE     , alpha_start = alpha,
                           auto_stop =  TRUE , stop_cond = stop_cond    , 
                           normalized_f = FALSE , ...);
    alpha_each[j] <- 0
    
    if (TRUE == result_PAFit$diverge_zero)
      alpha_each[j] <- -Inf
    else for (k in 1:length(cv_data$m_each)) 
      if (cv_data$m_each[k] != 0) {
        #PA              <- result_PAFit$A[cv_data$deg_each[k,] + 1] 
        #PA[is.na(PA)]   <- result_PAFit$A[length(result_PAFit$A)]
        PA              <-  cv_data$deg_each[k,]^alpha
        #print(PA)
        PA[PA == 0]     <- 1
        prob_PAFit      <- PA * result_PAFit$f[as.character(cv_data$stats$f_position)] 
        prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
        prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
        alpha_each[j] <- alpha_each[j] + 
                FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,as.character(cv_data$stats$f_position)] * 
                                                                cv_data$m_each[k])) 
      }
    #alpha_ea
    if (alpha_each[j] > max_val) {
      estimated_fitness <- result_PAFit$f[as.character(cv_data$stats$f_position)]
      s_optimal         <- rate_PAFit[j]
      max_val           <- alpha_each[j]
    }
  }
  #print(alpha_each)  
  
  if (print_out == TRUE)
    print(paste0("After first iteration: s_optimal = ",s_optimal))  
  

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
      result_PAFit <- PAFit(cv_data$stats, mode_f = "Log_linear", 
                            s = rate_PAFit[j], 
                            start_f      = estimated_fitness, 
                            alpha_start  = alpha, 
                            auto_stop    = TRUE,
                            only_f       = TRUE,
                            stop_cond    = stop_cond, 
                            normalized_f = FALSE,...);
      alpha_each[j] <- 0;
 
      if (TRUE == result_PAFit$diverge_zero)
        alpha_each[j] <- -Inf
      else for (k in 1:length(cv_data$m_each)) 
        if (cv_data$m_each[k] != 0) {
          PA              <-  cv_data$deg_each[k,]^alpha
          PA[PA == 0]     <- 1
          prob_PAFit      <- PA * result_PAFit$f[as.character(cv_data$stats$f_position)] 
          prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
          prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
          alpha_each[j]   <- alpha_each[j] + 
            FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,as.character(cv_data$stats$f_position)] * cv_data$m_each[k])) 
        }
      if (alpha_each[j] > max_val) {
        estimated_fitness <- result_PAFit$f[as.character(cv_data$stats$f_position)]
        s_optimal         <- rate_PAFit[j]
        max_val           <- alpha_each[j]
      }
    }
  if (print_out == TRUE)
    print(paste0("After second iteration: s_optimal = ",s_optimal))  
  
  
  rate_PAFit <- c(1/1.75 , 1/1.5 , 1/1.25 , 1.25 , 1.5 , 1.75) * s_optimal
  
  
  rate_PAFit <- rate_PAFit[rate_PAFit <= 10^4]
  
  rate_PAFit <- rate_PAFit[rate_PAFit >= 10^-1]
  
  alpha_each             <- rep(0,length(rate_PAFit))
  names(alpha_each)      <- rate_PAFit
  
  if (length(rate_PAFit) > 0)
    for (j in 1:length(rate_PAFit)) {
      #print(paste0("Processing case ",count, " of a maximum of ",total))
      result_PAFit <- PAFit(cv_data$stats, mode_f = "Log_linear", s = rate_PAFit[j], 
                            start_f = estimated_fitness, 
                            alpha_start = alpha,
                            only_f      = TRUE,
                            auto_stop   = TRUE, 
                            stop_cond = stop_cond, 
                            normalized_f = FALSE,...);
      alpha_each[j] <- 0;

      if (TRUE == result_PAFit$diverge_zero)
        alpha_each[j] <- -Inf
      else for (k in 1:length(cv_data$m_each)) 
          if (cv_data$m_each[k] != 0) {
              PA              <-  cv_data$deg_each[k,]^alpha
              PA[PA == 0]     <- 1
              prob_PAFit      <- PA * result_PAFit$f[as.character(cv_data$stats$f_position)] 
              prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
              prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
              alpha_each[j] <- alpha_each[j] + 
                   FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,as.character(cv_data$stats$f_position)] * cv_data$m_each[k])) 
        }
      
      if (alpha_each[j] > max_val) {
        estimated_fitness <- result_PAFit$f[as.character(cv_data$stats$f_position)]
        s_optimal         <- rate_PAFit[j]
        max_val           <- alpha_each[j]
      }
    }
  #print(alpha_each)  
  if (print_out == TRUE)
    print(paste0("Finally: s_optimal = ",s_optimal))  
  
  
  
  result    <- list(s_optimal         = s_optimal,
                    estimated_fitness = estimated_fitness)
  class(result) <- "CV_Result"
  
  return(result)
}


