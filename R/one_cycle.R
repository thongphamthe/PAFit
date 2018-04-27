.one_cycle <- function(cv_data                            ,
                       r                                  ,
                       s                                  ,
                       estimated_fitness_start = NULL     ,
                       estimated_PA_start      = NULL     ,
                       alpha_start             = NULL     ,
                       stop_cond               = 10^-6    ,
                       print_out               = FALSE    ,
                       cv_deg_thresh           = c(1,10)  ,
                       normal_start_f          = TRUE     ,
                       weight_f                = 0        ,
                       ...) { 
  
  FitMultinomial         <- function(true,dat){
    true[true == 0] <- 1
    return(sum(dat*log10(true)))
  }
 
  rate_PAFit             <- s
  
  if (is.null(estimated_PA_start)) {
      alpha_optimal        <- rep(-1,length(cv_deg_thresh))
  }
  
  alpha_each           <- matrix(0,nrow = length(cv_deg_thresh), ncol = length(rate_PAFit))
  colnames(alpha_each) <- rate_PAFit
  estimated_fitness    <- list(length = length(cv_deg_thresh))
  #estimated_PA        <- NULL
  max_val              <- rep(-Inf,length(cv_deg_thresh))
  s_optimal            <- rep(0,length(cv_deg_thresh))
  
  
  
  for (j in 1:length(rate_PAFit)) {
    
    #print(paste0("Processing case ",count, " of ",total))
    if (is.null(estimated_PA_start)) {
        if (!is.null(estimated_fitness_start)) { 
            result_PAFit  <- PAFit(cv_data$stats , mode_f    = "Log_linear"    , s = rate_PAFit[j], 
                                   auto_stop      =  TRUE , stop_cond  = stop_cond , normalized_f = FALSE ,
                                   normal_start_f = normal_start_f,
                                   weight_power   = weight_f, start_f = estimated_fitness_start,
                                   ...);
        }
        else 
            result_PAFit  <- PAFit(cv_data$stats     , mode_f    = "Log_linear"    , s = rate_PAFit[j], 
                                   auto_stop      =  TRUE , stop_cond  = stop_cond , normalized_f = FALSE ,
                                   normal_start_f = normal_start_f,
                                   weight_power   = weight_f, 
                                    ...);  
    }
    else {
        #print("inside here");  
        if (!is.null(estimated_fitness)) { 
            result_PAFit  <- PAFit(cv_data$stats     , true_A = estimated_PA_start, 
                                   s = rate_PAFit[j] , only_f = TRUE,
                                   auto_stop =  TRUE , 
                                   stop_cond  = stop_cond , normalized_f = FALSE ,
                                   normal_start_f = normal_start_f,
                                   weight_power   = weight_f, start_f = estimated_fitness_start,
                                   ...);  
            #plot(result_PAFit$center_k,result_PAFit$theta,log = "xy", main = "Inside loop");
            #print(paste0("Max of theta inside one cycle: ",max(estimated_PA_start)))
        }
        else 
            result_PAFit  <- PAFit(cv_data$stats     , true_A = estimated_PA_start, 
                                   s = rate_PAFit[j] , only_f = TRUE,
                                   auto_stop =  TRUE , 
                                   stop_cond  = stop_cond , normalized_f = FALSE ,
                                   normal_start_f = normal_start_f,
                                   weight_power   = weight_f,
                                   ...);    
    }
    if (is.null(estimated_PA_start)) {
        alpha_temp    <- result_PAFit$alpha
        #print(alpha_temp)
        if (TRUE == result_PAFit$diverge_zero) {
            alpha_each[,j] <- -Inf
        }
        else  { 
            for (dd in 1:length(cv_deg_thresh)) { 
                alpha_each[dd,j] <- 0  
                chosen_node_big  <- names(cv_data$stats$z_j[cv_data$stats$z_j >= cv_deg_thresh[dd]])
            for (k in 1:length(cv_data$m_each)) 
                if (cv_data$m_each[k] != 0) {
                #PA              <- result_PAFit$A[cv_data$deg_each[k,] + 1] 
                #PA[is.na(PA)]   <- result_PAFit$A[length(result_PAFit$A)]
                chosen_node     <- chosen_node_big[cv_data$deg_each[k,chosen_node_big] != 0]   
                PA              <-  cv_data$deg_each[k,chosen_node]^alpha_temp
                #print(PA)
                if (sum(PA == 0) > 0)
                    print("This should not happen")  
                #PA[PA == 0]     <- 1
                fitness         <- rep(1,dim(cv_data$deg_each[,chosen_node, drop = FALSE])[2])
                names(fitness)  <- colnames(cv_data$deg_each[,chosen_node, drop = FALSE])
                fitness[chosen_node] <- result_PAFit$f[chosen_node] 
            
                prob_PAFit      <- PA * fitness
                prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
                prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
                alpha_each[dd,j]   <- alpha_each[dd, j] + 
                                      FitMultinomial(true = as.vector(prob_PAFit), 
                                                     dat = as.vector(cv_data$prob_em_each[k,chosen_node] * 
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
    } else {
         ### estimated_PA is provided
          if (TRUE == result_PAFit$diverge_zero) {
              alpha_each[,j] <- -Inf
          } else  { 
             
              for (dd in 1:length(cv_deg_thresh)) { 
                  alpha_each[dd,j] <- 0  
                  chosen_node_big      <- names(cv_data$stats$z_j[cv_data$stats$z_j >= cv_deg_thresh[dd]])
                  for (k in 1:length(cv_data$m_each)) 
                  if (cv_data$m_each[k] != 0) {
                      chosen_node      <- chosen_node_big[!is.na(result_PAFit$A[cv_data$deg_each[k,chosen_node_big] + 1])]
                      chosen_node      <- chosen_node[result_PAFit$A[cv_data$deg_each[k,chosen_node] + 1] != 0]
                      
                      #PA              <- result_PAFit$A[cv_data$deg_each[k,chosen_node] + 1] 
                      #PA              <- result_PAFit$A[cv_data$deg_each[k,chosen_node] + 1]
                      PA              <- (cv_data$deg_each[k,chosen_node])^alpha_start
                      if (sum(is.na(PA)) > 0) {
                          #PA[is.na(PA)] <- max(estimated_PA_start)  
                          #PA[is.na(PA)]   <- cv_data$deg_each[k,chosen_node][is.na(PA)]^alpha_start
                          print("It should not happen")
                      }
                      if (sum(PA == 0) > 0)
                          PA[PA == 0]     <- 1
                      fitness         <- rep(1,dim(cv_data$deg_each[,chosen_node, drop = FALSE])[2])
                      names(fitness)  <- colnames(cv_data$deg_each[,chosen_node, drop = FALSE])
                      fitness[chosen_node] <- result_PAFit$f[chosen_node] 
              
                      prob_PAFit      <- PA * fitness
                      prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
                      prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
                      alpha_each[dd,j] <- alpha_each[dd, j] + 
                                             FitMultinomial(true = as.vector(prob_PAFit), 
                                                           dat = as.vector(cv_data$prob_em_each[k,chosen_node] * 
                                                           cv_data$m_each[k])) 
                 }
            }
          }
          for (dd in 1:length(cv_deg_thresh)) 
              if (alpha_each[dd,j] > max_val[dd]) {
                  estimated_fitness[[dd]]  <- result_PAFit$f[as.character(cv_data$stats$f_position)]
                  s_optimal[dd]            <- rate_PAFit[j]
                  max_val[dd]              <- alpha_each[dd,j]
                  #alpha_optimal[dd]        <- result_PAFit$alpha
              }
          }
  }
  
  
  
  ok_index                <- 1:length(cv_deg_thresh)
  
  s_optimal_final         <- mean(s_optimal[ok_index])
  
  #s_optimal_final         <- 10^(1/length(s_optimal[ok_index])* log10(prod(s_optimal[ok_index])))
  
  if (is.null(estimated_PA_start))
      alpha_optimal_final     <- mean(alpha_optimal[ok_index])
  estimated_fitness_final <- estimated_fitness[[ok_index[1]]]
  
  if (length(ok_index) > 1)
    for (temp_jj in 2:length(ok_index))
      estimated_fitness_final <- estimated_fitness_final + estimated_fitness[[ok_index[temp_jj]]] 
  estimated_fitness_final <- estimated_fitness_final / length(ok_index)
  
  

  ### One pass to find optimal r #######

  chosen_node_big <- names(cv_data$stats$z_j[cv_data$stats$z_j >= cv_deg_thresh[1]])

  ratio_vec_PAFit        <- sort(r,decreasing = TRUE)
  PA_each                <- rep(0,length(ratio_vec_PAFit))
  names(PA_each)         <- ratio_vec_PAFit
  max_val                <- -Inf 
  
  
  estimated_PA           <- NULL
  
  switch_flag            <- 0
  
  for (i in 1:length(ratio_vec_PAFit)) {
   
    #print(paste0("Processing case ",count, " of a maximum of ",total))
    if (!is.null(estimated_PA))
      result_PAFit <- PAFit(cv_data$stats, 
                            s = s_optimal_final, 
                            r = ratio_vec_PAFit[i], 
                            only_PA   = TRUE, 
                            true_f    = estimated_fitness_final, 
                            auto_stop = TRUE, 
                            start_A   = estimated_PA,
                            stop_cond = stop_cond, 
                            normal_start_f = normal_start_f,
                            weight_power   = weight_f ,
                            normalized_f = FALSE, ...)
    
    else {
        if (is.null(estimated_PA_start)) {
            result_PAFit <- PAFit(cv_data$stats, 
                                  s = s_optimal_final, 
                                  r = ratio_vec_PAFit[i], 
                                  only_PA = TRUE, 
                                  true_f = estimated_fitness_final, 
                                  auto_stop =  TRUE, 
                                  alpha_start = alpha_optimal_final,
                                  stop_cond = stop_cond, 
                                  normal_start_f = normal_start_f,
                                  normalized_f = FALSE, ...)
        }
        else 
          result_PAFit <- PAFit(cv_data$stats, 
                                s = s_optimal_final, 
                                r = ratio_vec_PAFit[i], 
                                only_PA = TRUE, 
                                true_f = estimated_fitness_final, 
                                auto_stop =  TRUE, 
                                start_A = estimated_PA_start,
                                stop_cond = stop_cond, 
                                normal_start_f = normal_start_f,
                                normalized_f = FALSE, ...)  
    }
    
    
    #alpha_temp    <- result_PAFit$alpha  
    
    
    if (TRUE == result_PAFit$diverge_zero)
        PA_each[i] <- -Inf
    else for (k in 1:length(cv_data$m_each))
        if (cv_data$m_each[k] != 0) { 
            chosen_node    <- chosen_node_big[!is.na(result_PAFit$A[cv_data$deg_each[k,chosen_node_big] + 1])] 
            chosen_node    <- chosen_node[result_PAFit$A[cv_data$deg_each[k,chosen_node] + 1] != 0]
            #PA              <- cv_data$deg_each[k,chosen_node]^alpha_temp
            PA              <- result_PAFit$A[cv_data$deg_each[k,chosen_node] + 1]
            if (sum(is.na(PA)) > 0) {
                print("it should not happen here")  
                PA[is.na(PA)]   <- cv_data$deg_each[k,chosen_node][is.na(PA)]^alpha_temp
            }
            #print(PA)
            #PA[PA == 0]     <- mean(result_PAFit$A)
            if (sum(PA == 0) > 0) {
                print("This should not happen")  
                PA[PA == 0]     <- cv_data$deg_each[k,chosen_node][PA == 0]^alpha_temp
            }
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
         alpha_op_final    <- result_PAFit$alpha
    } else if (PA_each[i] > max_val) {
          lambda_optimal    <- result_PAFit$lambda  
          r_optimal         <- ratio_vec_PAFit[i]
          max_val           <- PA_each[i]
          estimated_PA      <- result_PAFit$theta  
          alpha_op_final    <- result_PAFit$alpha
    }
  }
  #if (print_out == TRUE)
  
  #print(paste0("s_optimal = ",s_optimal_final, " ,r_optimal = ", r_optimal)) 
  
  result    <- list(r_optimal         = r_optimal, 
                    lambda_optimal    = lambda_optimal, 
                    s_optimal         = s_optimal_final,
                    estimated_fitness = estimated_fitness_final,
                    estimated_PA      = estimated_PA,
                    cv_deg_thresh     = cv_deg_thresh,
                    alpha_optimal     = alpha_op_final,
                    
                    center_k          = result_PAFit$center_k)
  class(result) <- "CV_Result"
  return(result)
}