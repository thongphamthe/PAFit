# function to estimate jointly the attachment function and node fitness  
PAFit <- function(net_stat, 
                  only_PA        = FALSE       , only_f         = FALSE        , 
                  mode_f         = "Linear_PA" ,
                  true_A         = NULL        , true_f         = NULL         , 
                  
                  mode_reg_A     = 0           , weight_PA_mode = 1            ,
                  s              = 10          , lambda         = 1            , 
                  auto_lambda    = TRUE        , r              = 0.01         ,
                  
                  alpha_start    = 0.5         , start_mode_A   = "Log_linear" , 
                  start_mode_f   = "Constant"  , start_A        = NULL         ,
                  start_f        = NULL        ,
                  
                  auto_stop      = TRUE        , stop_cond      = 10^-8        , 
                  iteration      = 200         , max_iter       = 200000       , 
                  debug          = FALSE       , q              = 1            , 
                  step_size      = 0.5         ,
                  
                  normalized_f   = FALSE       , interpolate    = FALSE        ,
                  weight_power   = 1           , normal_start_f = FALSE        ) {
  if ((net_stat$only_PA == TRUE) & (only_PA == FALSE)) {
    warning("The net_stat does not support estimation of node fitness. It will be run with option 'only_PA = TRUE'.
            Please re-run GetStatistics again with the option 'only_PA = FALSE' if you also want to estimate fitnesses.")
    only_PA <- TRUE  
  }
  
  

  oopts <- options(scipen = 999)
  on.exit(options(oopts))
  shape          <- s
  rate           <- s
  ratio          <- r
  if (s <= 0)
    stop("Error: shape and rate should be positive number")
  
  non_zero_theta     <- which(net_stat$sum_m_k > 0)
  num_nonzero        <- length(non_zero_theta)
  if ((num_nonzero == 1) && (only_f == FALSE)) {
     # only one non-zero bin
      stop(paste0("Error: There is only one bin that has a non-zero number of new edges (bin ",which(net_stat$sum_m_k > 0),"). To estimate the PA function, we need at least two bins with non-zero number of new edges."))  
  }
  if ((num_nonzero == 0) && (only_f == FALSE)) {
    # no non-zero bin
    stop(paste0("Error: There is no bin that has a non-zero number of new edges. To estimate the PA function, we need at least two bins with non-zero number of new edges."))  
  }
  theta              <- rep(0,length(net_stat$sum_m_k))
  
  
  for (ii in 1:length(theta))
    if ((mode_f[1] == "Constant_PA") & (only_f == TRUE)) {
      theta[ii] <- 1;  
    } else if ((only_f == TRUE) & (!is.null(true_A))) {
      
      theta[ii] <- true_A[ii];    
    }
    else if (!is.null(start_A)) {
        theta[ii] <- start_A[ii];  
     } else if ("Log_linear" != mode_f) {
    if ("Log_linear" == start_mode_A[1]){    
      if (ii <= net_stat$start_deg & ii %in% non_zero_theta) {
        theta[ii] <- ii^alpha_start
      }
      else 
        if (ii > net_stat$start_deg & ii %in% non_zero_theta) {
         
          if (net_stat$begin_deg[ii - net_stat$start_deg] > 0)  
            theta[ii] <- net_stat$begin_deg[ii - net_stat$start_deg]^alpha_start
          else     theta[ii] <- 1 
        }
    } else if ("Random" == start_mode_A[1]) {
      theta <- runif(length(theta),0,1);  
    }
  } else if ("Log_linear" == mode_f) {
    if (ii <= net_stat$start_deg & ii %in% non_zero_theta) {
      theta[ii] <- ii
    }
    else 
      if (ii > net_stat$start_deg & ii %in% non_zero_theta) {
      
        if (net_stat$begin_deg[ii - net_stat$start_deg] > 0)  
          theta[ii] <- net_stat$begin_deg[ii - net_stat$start_deg]
        else     theta[ii] <- 1 
      }
  }
 

  alpha <- 1.0
  if ("Log_linear" != mode_f)
      if (length(theta) > 1)  {
          if (theta[length(theta)] == 1 || theta[length(theta)] == 0)
              theta[length(theta)] <- theta[length(theta) - 1] + 1
  }
 
  theta[which(theta <= 0)] <- 1
  if (mode_f[1] != "Log_linear")
    theta[-non_zero_theta] <- 0
  if (mode_f[1] == "Log_linear") {
    alpha <- 1.0 * alpha_start; #starting value of alpha
  }
  #else theta         <- theta/sum(theta)
 
  #if (include_zero == 0)
  non_zero_f    <- which(net_stat$z_j > 0)
  
  #if (shape <= 1 & rate <= 1)
  #    non_zero_f <- which(net_stat$z_j > 0)
  
  if (!is.null(true_f)) {
      f <- true_f  
  }
  else if (TRUE == only_PA) {
      f <- rep(1,length(net_stat$f_position)) 
  } else if (!is.null(start_f)) {
      f <- start_f  
  } else if ("Constant" == start_mode_f) {
      
    
      f <- rep(1,length(net_stat$f_position)) 
      
      # reasoning: the final nodes with z_j = 0 should have fitness 1
      # earlier nodes with z_j = 0 should have fitness less than 1
  }
  else if ("Random" == start_mode_f) {
      f <- rgamma(n = length(net_stat$f_position), shape = s, rate = s);
      f <- length(f) * f / sum(f)
  }
  
  
  
  f[net_stat$z_j == 0] <- 1
  
  if (TRUE == only_PA) {
    if (is.null(true_f))
      f[] <- 1
    else
      f[] <- true_f  
  }

  
  max_time <- max(net_stat$appear_time)
  weight_f <- max_time/ (max_time - net_stat$appear_time[as.character(as.integer(net_stat$f_position))] + 1)
 
  weight_f <- weight_f^0
  #weight_f <- rep(1,length(net_stat$f_position))
  
 
  
  
  log10_likelihood   <- vector()
  # get the center of each bin
  center_k     <- net_stat$center_k
  
  
  if (mode_reg_A != 0)
      mode_reg_A <- 2  
  
  
  if (0 == mode_reg_A) {
    update_theta   <- non_zero_theta[-c(1 , 2 , num_nonzero , num_nonzero - 1)]
    noupdate_theta <- non_zero_theta[c(1 , 2 , num_nonzero , num_nonzero - 1)]
    if (num_nonzero > 2)
      minus_1        <- non_zero_theta[-c(1 , num_nonzero , num_nonzero - 1 , num_nonzero - 2)]
    if (num_nonzero > 3)
      minus_2    <- non_zero_theta[-c(num_nonzero , num_nonzero - 1 , num_nonzero - 2 , num_nonzero - 3)] 
    
    plus_1         <- non_zero_theta[-c(1 , 2 , 3 , num_nonzero)]
    
    if (length(non_zero_theta) > 4)
      plus_2         <- non_zero_theta[-c(1 , 2 , 3 , 4)] 
    #num_ok_log10     <- 1
  } else if (1 == mode_reg_A) {
    ok_log10     <- which(center_k[non_zero_theta] > 1)
    not_ok_log10 <- which(center_k[non_zero_theta] <= 1)
    num_ok_log10 <- length(ok_log10) 
    
    update_theta   <- non_zero_theta[ok_log10[-c(1, 2, num_ok_log10, num_ok_log10 - 1)]]
    noupdate_theta <- non_zero_theta[ok_log10[c(1, 2, num_ok_log10, num_ok_log10 - 1)]]
    if (num_ok_log10 > 2)
      minus_1 <- non_zero_theta[ok_log10[-c(1, num_ok_log10, num_ok_log10 - 1, num_ok_log10 - 2)]]
    if (num_ok_log10 > 3)
      minus_2 <- non_zero_theta[ok_log10[-c(num_ok_log10, num_ok_log10 - 1, num_ok_log10 - 2, num_ok_log10 - 3)]] 
    
    plus_1  <- non_zero_theta[ok_log10[-c(1, 2, 3, num_ok_log10)]]
    
    if (length(ok_log10) > 4)
      plus_2  <- non_zero_theta[ok_log10[-c(1, 2, 3, 4)]] 
    
  } else if (2 == mode_reg_A) {
    ok_log10     <- which(center_k[non_zero_theta] > 0)
    not_ok_log10 <- which(center_k[non_zero_theta] == 0)
    num_ok_log10 <- length(ok_log10) 

    if (num_ok_log10 >= 3) { 
      update_theta   <- non_zero_theta[ok_log10[-c(1, 2, num_ok_log10, num_ok_log10 - 1)]]
      noupdate_theta <- non_zero_theta[ok_log10[c(1, 2, num_ok_log10, num_ok_log10 - 1)]]
      minus_1 <- non_zero_theta[ok_log10[-c(1, num_ok_log10, num_ok_log10 - 1, num_ok_log10 - 2)]]
      minus_2 <- non_zero_theta[ok_log10[-c(num_ok_log10, num_ok_log10 - 1, num_ok_log10 - 2, num_ok_log10 - 3)]] 
      plus_1  <- non_zero_theta[ok_log10[-c(1, 2, 3, num_ok_log10)]]
      plus_2  <- non_zero_theta[ok_log10[-c(1, 2, 3, 4)]] 
      u       <- matrix(0,nrow = num_ok_log10, ncol = num_ok_log10)

      
      for (k1 in 2:(num_ok_log10 - 1)) {
        u[k1,k1]       <- 1/(log(center_k[non_zero_theta[ok_log10[k1 + 1 ]]]) - log(center_k[non_zero_theta[ok_log10[k1]]])) +
          1/(log(center_k[non_zero_theta[ok_log10[k1]]]) - log(center_k[non_zero_theta[ok_log10[k1 - 1]]])) 
        u[k1,k1 - 1]   <- 1/(log(center_k[non_zero_theta[ok_log10[k1]]]) - log(center_k[non_zero_theta[ok_log10[k1 - 1]]]))
        u[k1,k1 + 1]   <- 1/(log(center_k[non_zero_theta[ok_log10[k1 + 1 ]]]) - log(center_k[non_zero_theta[ok_log10[k1]]]))
      }
      
      extract_u <- function(i,j) {
        result <- NULL
      
        for (index in 1:length(i))
          result <- c(result,u[i[index],j[index]])  
        return(result)
      }
      update_k <- 3:(num_ok_log10-2)
    }
    
  }
  
  PA_offset <- 1
  #if(theta[1] != 0)
  theta[1]  <- PA_offset
  
  #starting value for the offset
  offset <- 1
  #weights of the regularization term
  if (weight_PA_mode[1] == 0)
    w_k <- net_stat$sum_m_k/sum(net_stat$sum_m_k)
  else {
    if (weight_PA_mode[1] == 1)
      w_k <- rep(1/length(theta),length(theta))
    else if (weight_PA_mode[1] == 2)
      w_k <- 1/net_stat$sum_m_k/sum(1/net_stat$sum_m_k)
  }
  
  if (TRUE == auto_lambda) {
    #lambda <- ratio   
    lambda <- ratio * sum(net_stat$sum_m_k)
  }
  if (TRUE == auto_stop)
    iteration <- max_iter
  
  normalized_const <- rep(0, dim(net_stat$n_tk)[1])
  alpha_series <- c()
  
  break_flag  <- FALSE
  count_break <- 0
  diverge_zero         <- FALSE
  diverge_zero_theta   <- FALSE
  non_zero_theta_start <- theta > 10^-30
  non_zero_f_start     <- f     > 10^-30
  
  cal_reg_A <- function() {
    
    if (lambda > 0) {  
      if (0 == mode_reg_A) {
        return(sum(lambda * w_k[non_zero_theta[-c(1,num_nonzero)]] * (log(theta[non_zero_theta[-c(1 , 2)]]) + 
                                                                        log(theta[non_zero_theta[-c(num_nonzero , num_nonzero - 1)]]) - 
                                                                        2 * log(theta[non_zero_theta[-c(1 , num_nonzero)]]))^2))
      } else if (1 == mode_reg_A) {
        return(sum(lambda * w_k[non_zero_theta[ok_log10[-c(1 , num_ok_log10)]]] *
                     (log(theta[non_zero_theta[ok_log10[-c(1 , 2)]]]) / log(center_k[non_zero_theta[ok_log10[-c(1 , 2)]]]) + 
                        log(theta[non_zero_theta[ok_log10[-c(num_ok_log10 , num_ok_log10 - 1)]]]) / 
                        log(center_k[non_zero_theta[ok_log10[-c(num_ok_log10 , num_ok_log10 - 1)]]])  - 
                        2 * log(theta[non_zero_theta[ok_log10[-c(1 , num_ok_log10)]]]) / 
                        log(center_k[non_zero_theta[ok_log10[-c(1 , num_ok_log10)]]]) ) ^ 2))
      } else if (2 == mode_reg_A) {
       
        if (num_ok_log10 >= 3) {
          return_temp <- sum(lambda * w_k[non_zero_theta[ok_log10[-c(1,num_ok_log10)]]] *
                               ( (
                                 (log(theta[non_zero_theta[ok_log10[-c(1,2)]]]) -  log(theta[non_zero_theta[ok_log10[-c(1,num_ok_log10)]]])) / 
                                   (log(center_k[non_zero_theta[ok_log10[-c(1,2)]]]) -  log(center_k[non_zero_theta[ok_log10[-c(1,num_ok_log10)]]]))
                               ) -  
                                 (
                                   (log(theta[non_zero_theta[ok_log10[-c(1,num_ok_log10)]]]) -  log(theta[non_zero_theta[ok_log10[-c(num_ok_log10,num_ok_log10 - 1)]]])) / 
                                     (log(center_k[non_zero_theta[ok_log10[-c(1,num_ok_log10)]]]) -  log(center_k[non_zero_theta[ok_log10[-c(num_ok_log10, num_ok_log10 - 1)]]])) 
                                 ) 
                               )^2)
          if (length(not_ok_log10) > 0)
            return_temp <- return_temp +  lambda * w_k[non_zero_theta[not_ok_log10]] * 
              (log(theta[non_zero_theta[not_ok_log10]]) - log(theta[non_zero_theta[ok_log10[1]]])) ^2    
          return(return_temp);  
        } else return(0);
        
      }
    } else return(0);
  }
  
  objective_function_value <- function(theta,f,offset,normalized_const) {
    time_non_zero       <- which(normalized_const > 0)
    non_zero_f_temp     <- which(f > 0)
    non_zero_theta_temp <- which(theta > 0) 
  
    if (TRUE == only_f) {
      value <- sum(net_stat$z_j[non_zero_f_temp] * log(f[non_zero_f_temp])) - 
        sum(net_stat$m_t[time_non_zero] * log(normalized_const[time_non_zero])) + 
        sum(net_stat$offset_m_tk) * log(offset) +   
        (sum((shape / weight_f[non_zero_f_temp] - 1) * log(f[non_zero_f_temp])) - 
           sum(rate / weight_f[non_zero_f_temp]  * (f[non_zero_f_temp])));  
    } else if (FALSE == only_PA) {
      value <- sum(net_stat$z_j[non_zero_f_temp] * log(f[non_zero_f_temp])) + 
        sum(net_stat$offset_m_tk)*log(offset) +  
        (shape - 1) * log(offset) - rate * offset +    
        sum(net_stat$sum_m_k[non_zero_theta_temp] * log(theta[non_zero_theta_temp])) - 
        sum(net_stat$m_t[time_non_zero] * log(normalized_const[time_non_zero])) + 
        (sum((shape / weight_f[non_zero_f_temp] - 1) * log(f[non_zero_f_temp])) - 
           sum(rate / weight_f[non_zero_f_temp]  * (f[non_zero_f_temp]))) - cal_reg_A()
     
    }  else if ((TRUE == only_PA) && (!is.null(true_f))) {
      value <- sum(net_stat$z_j[non_zero_f_temp] * log(f[non_zero_f_temp])) + 
        sum(net_stat$sum_m_k[non_zero_theta_temp] * log(theta[non_zero_theta_temp])) - 
        sum(net_stat$m_t[time_non_zero] * log(normalized_const[time_non_zero])) - cal_reg_A()
      
    }  else if (is.null(true_f) && (TRUE == only_PA))  { 
     
      value <- sum(net_stat$sum_m_k[non_zero_theta_temp] * log(theta[non_zero_theta_temp])) - 
        sum(net_stat$m_t[time_non_zero] * log(normalized_const[time_non_zero]))
      
      
      value <- value - cal_reg_A()
    }
    
    return(value)
  }
  
  if (q > 1) {
    parameter_save <- matrix(0,nrow = length(theta) + length(f) + 1, ncol = q)
    U              <- matrix(0,nrow = length(theta) + length(f) + 1, ncol = q - 1)
    V              <- matrix(0,nrow = length(theta) + length(f) + 1, ncol = q - 1)
    candidate      <- rep(0,length(theta) + length(f) + 1)              # theta, f, and offset values
  }
  
  candidate_ok     <- 0
  candidate_accept <- 0
  
  
  for (i in 1:iteration) {
    
    if (mode_f[1] != "Log_linear") {
      if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))) {
        .normalized_constant(normalized_const,net_stat$node_degree,theta,f,net_stat$offset_tk,offset) 
      } else if (is.null(true_f) && (TRUE == only_PA)) {
        normalized_const <- as.vector(net_stat$n_tk[,non_zero_theta]%*%theta[non_zero_theta])
      }
   
      
      log10_likelihood        <- c(log10_likelihood,objective_function_value(theta,f,offset,normalized_const))  ;
      names(log10_likelihood) <- NULL
      time_non_zero           <- which(normalized_const != 0)
      non_zero_f_temp         <- which(f > 0)
      
      
      if (TRUE == debug) {
        print(log10_likelihood[length(log10_likelihood)])
      }
      
      break_flag <- FALSE;
      if (length(log10_likelihood) > 1) {
        if (abs(log10_likelihood[length(log10_likelihood)] - log10_likelihood[length(log10_likelihood) - 1]) /
            abs(log10_likelihood[length(log10_likelihood)]) > 10^-15)  
          if (log10_likelihood[length(log10_likelihood)] < log10_likelihood[length(log10_likelihood) - 1]) {
                       #stop("Warning: log likelihood decreased.")  
            break_flag <- TRUE;
          }
      }
      
      if (TRUE == auto_stop)
        if (length(log10_likelihood) > 1)
          tryCatch({if (abs(log10_likelihood[length(log10_likelihood)] - log10_likelihood[length(log10_likelihood) - 1]) / 
                        (abs(log10_likelihood[length(log10_likelihood) - 1]) + 1) < stop_cond)
            break_flag <- TRUE;},error = function(e) { 
              break_flag <- TRUE;})
      
      
      ######################### quasi-Newton acceleration #########################
      ####  calculate the smallest approximation to the second derivative at current point ####
      if ((FALSE == only_f) && (FALSE == only_PA)) {
        if (q > 1) {
          flag <- 0  
          if (i > q) {
            current_pos <- c(theta,f,offset)
            U <- parameter_save[,1:(q-1)] - parameter_save[,q]
            V <- parameter_save[,2:q]     - current_pos            
            
            if ((FALSE == only_f) && (FALSE == only_PA)) {
              candidate <- tryCatch(current_pos + step_size * ifelse(q > 2, V%*%solve(crossprod(U,U) - 
                                                                                        crossprod(U,V),crossprod(U, V[,q-1])), 1/(sum(U* U) - sum(U*V)) * sum(V*U) * V), 
                                    error = function(e) {
                                      flag <- 1; return(current_pos)})
              positive <- prod(candidate > 0)      
            }
            
            if ((0 == flag) && (1 == positive)) {
              candidate_ok <- candidate_ok + 1
              
              theta_temp                 <- candidate[1:length(theta)]
              f_temp                     <- candidate[(length(theta) + 1):(length(theta) + length(f))]
              offset_temp                <- candidate[length(candidate)]  
              normalized_const_temp      <- rep(0, dim(net_stat$n_tk)[1])
              
              if (FALSE == only_PA) {
                .normalized_constant(normalized_const_temp,net_stat$node_degree,theta_temp,f_temp,net_stat$offset_tk,offset_temp)
              }  
              else if (TRUE == only_PA)
                normalized_const_temp <- as.vector(net_stat$n_tk[,non_zero_theta]%*% theta_temp[non_zero_theta])
              
              time_temp <- which(normalized_const_temp != 0)   
              log10_candidate <- objective_function_value(theta_temp,f_temp,offset_temp,normalized_const_temp)
              if (log10_candidate > log10_likelihood[length(log10_likelihood)]) {
                  theta  <- candidate[1 : length(theta)]
                  f      <- candidate[(length(theta) + 1) : (length(theta) + length(f))]
                  offset <- candidate[length(candidate)]  
                  log10_likelihood[length(log10_likelihood)]  <- log10_candidate
                  normalized_const                        <- normalized_const_temp
              }
              else {
               
              }
            }
            parameter_save[,1:(q-1)]             <- parameter_save[,2:q] 
            parameter_save[1:length(theta),q]    <- theta
            parameter_save[(length(theta) + 1):(length(theta) + length(f)),q] <- f
            parameter_save[dim(parameter_save)[1],q] <- offset
          }
          else {
            parameter_save[1:length(theta),i]    <- theta
            parameter_save[(length(theta) + 1):(length(theta) + length(f)),i] <- f
            parameter_save[dim(parameter_save)[1],i] <- offset
          }

          if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))){
            .normalized_constant(normalized_const,net_stat$node_degree,theta,f,net_stat$offset_tk,offset) 
          }
          else normalized_const <- as.vector(net_stat$n_tk[,non_zero_theta]%*%theta[non_zero_theta])
          #.normalized_constant(normalized_const, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
          
          
        }
      }
      ############### End of quasi-Newton acceleration ################
      
      
      if (break_flag) {
        break;   
      } else { 
        break_flag  <- FALSE
      }
      
      
      
      
      
      ##################### Update f ######################
      if ((FALSE == only_PA) && (is.null(true_f))) {
        .update_f_new(f,non_zero_f,net_stat$node_degree,theta,net_stat$z_j,normalized_const,net_stat$m_t,
                  shape,rate,offset,weight_f)
        
      
        
        ### NORMALIZED f
        #normalize_f <- mean(c(f,offset))
        #f           <- f/normalize_f
        #offset      <- offset/normalize_f
        
        non_zero_f_now <- f > 10^-30
        
        if (length(non_zero_f_now) > 0)
          for (lll in 1:length(non_zero_f_now))
            if (non_zero_f_start[lll] != non_zero_f_now[lll]) {
              
              diverge_zero <- TRUE
              break; 
            }
        if (TRUE == diverge_zero)
          break;   
        
        
        #update offset
        #if (net_stat$deg_thresh > 0)
        if (FALSE == only_f) {
          .normalized_constant(normalized_const, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
          time_non_zero     <- which(normalized_const != 0)
        }
        #offset <- .update_offset(net_stat$offset_tk, net_stat$offset_m_tk, theta, normalized_const,net_stat$m_t, shape,rate)
     
        #.normalized_constant(normalized_const, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
        #
       
      }
      
      #### Update PA_offset if mode_f = "Linear_PA" #####
      #if ((TRUE == only_f) && (mode_f[1] == "Linear_PA")) {
      #    PA_offset <- .update_PA_offset(normalized_const,f,net_stat$node_degree,net_stat$m_t,net_stat$sum_m_k,
      #                                   net_stat$offset_tk);
      #    theta[1]  <- PA_offset
      #}
      #####################  Update A #######################
      
      if (FALSE == only_f) {
        if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))) {
          
          temp5  <- .coeff_theta(net_stat$node_degree, f, normalized_const,net_stat$m_t,net_stat$start_deg + net_stat$g)   
          if (length(time_non_zero) > 1 && length(non_zero_theta) > 1) {
            temp4  <- temp5[non_zero_theta] + colSums(net_stat$m_t[time_non_zero] / normalized_const[time_non_zero] * offset *
                                                        net_stat$offset_tk[time_non_zero,non_zero_theta, drop = FALSE])
          } else {
            temp4  <- temp5[non_zero_theta] + net_stat$m_t[time_non_zero] / normalized_const[time_non_zero] * offset *
              net_stat$offset_tk[time_non_zero,non_zero_theta]  
          }
        }
        else {
          if (length(time_non_zero) > 1 && length(non_zero_theta) > 1) {  
             
              temp4 <- colSums(net_stat$m_t[time_non_zero]/normalized_const[time_non_zero] * 
                               net_stat$n_tk[time_non_zero,non_zero_theta, drop = FALSE])
          } else {
            temp4 <- net_stat$m_t[time_non_zero]/normalized_const[time_non_zero] * 
                       net_stat$n_tk[time_non_zero,non_zero_theta]  
          }
        }
       
        if ((lambda <= 0) || ((0 != mode_reg_A) && (num_ok_log10 < 3))) {
          
          theta[non_zero_theta] <- net_stat$sum_m_k[non_zero_theta] / temp4
        }
        else if (0 == mode_reg_A) {
          g_1  <- function(x) {
            (net_stat$sum_m_k[non_zero_theta[1]] - 2*w_k[non_zero_theta[2]] * lambda * (-2 * log(theta[non_zero_theta[2]]) + 
                                                                                          log(theta[non_zero_theta[3]]) - 3 * log(theta[non_zero_theta[1]])))/x -   
              temp4[1] - 8 * lambda * w_k[non_zero_theta[2]] * log(x) / x}
          if (length(non_zero_theta) > 3) {
            g_2 <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[2]] - 2 * lambda * 
                 (2 * w_k[non_zero_theta[2]] * (-2 * log(theta[non_zero_theta[2]]) -log(theta[non_zero_theta[3]]) - log(theta[non_zero_theta[1]])) + 
                    w_k[non_zero_theta[3]]*(-2*log(theta[non_zero_theta[3]]) +log(theta[non_zero_theta[4]]) -3 * log(theta[non_zero_theta[2]]) ))) / x - 
                temp4[2] - (16 * w_k[non_zero_theta[2]] + 8 * w_k[non_zero_theta[3]]) * lambda * log(x) / x}
          } else if (length(non_zero_theta) == 3) {
            g_2 <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[2]] - 2 * lambda * 
                 (2*w_k[non_zero_theta[2]] * (-2 * log(theta[non_zero_theta[2]]) -log(theta[non_zero_theta[3]]) - log(theta[non_zero_theta[1]])) + 
                    w_k[non_zero_theta[3]]*(-2*log(theta[non_zero_theta[3]]) -  3 * log(theta[non_zero_theta[2]]) ))) / x - 
                temp4[2] - (16 * w_k[non_zero_theta[2]] + 8 * w_k[non_zero_theta[3]]) * lambda * log(x) / x}          
          }
          if (num_nonzero > 3) {  
            g_semiend <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[num_nonzero - 1]] - 2*lambda * 
                 (2*w_k[non_zero_theta[num_nonzero - 1]]*(-2*log(theta[non_zero_theta[num_nonzero - 1]]) - log(theta[non_zero_theta[num_nonzero]]) - log(theta[non_zero_theta[num_nonzero - 2]])) +
                    w_k[non_zero_theta[num_nonzero -2]] *(-3*log(theta[non_zero_theta[num_nonzero - 1]]) + log(theta[non_zero_theta[num_nonzero - 3]]) - 2*log(theta[non_zero_theta[num_nonzero - 2]])))) / x - 
                temp4[num_nonzero - 1] - 
                (16 * w_k[non_zero_theta[num_nonzero - 1]] + 8*w_k[non_zero_theta[num_nonzero -2]]) * lambda * log(x) / x}
          } else if (num_nonzero == 3) {
            g_semiend <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[num_nonzero - 1]] - 2 * lambda * 
                 (2 * w_k[non_zero_theta[num_nonzero - 1]] * (-2 * log(theta[non_zero_theta[num_nonzero - 1]]) - log(theta[non_zero_theta[num_nonzero]]) - log(theta[non_zero_theta[num_nonzero - 2]])) +
                    w_k[non_zero_theta[num_nonzero -2]] *(-3 * log(theta[non_zero_theta[num_nonzero - 1]]) - 2 * log(theta[non_zero_theta[num_nonzero - 2]])))) / x - 
                temp4[num_nonzero - 1] - 
                (16 * w_k[non_zero_theta[num_nonzero - 1]] + 8 * w_k[non_zero_theta[num_nonzero -2]]) * lambda * log(x) / x}   
            
          }
          
          g_end     <-  function(x) {
            (net_stat$sum_m_k[non_zero_theta[num_nonzero]] - 2*lambda * w_k[non_zero_theta[num_nonzero -1]] *(-3*log(theta[non_zero_theta[num_nonzero]]) + log(theta[non_zero_theta[num_nonzero - 2]]) - 2*log(theta[non_zero_theta[num_nonzero - 1]]))) / x  -   
              temp4[num_nonzero] - 8*w_k[non_zero_theta[num_nonzero -1]]*lambda*log(x)/x}  
          
          if (length(update_theta) > 0) {
            U_1 <- net_stat$sum_m_k[update_theta] - 2 * lambda * (2*w_k[update_theta]*(-2*log(theta[update_theta]) - log(theta[plus_1]) - log(theta[minus_1])) + 
                                                                    w_k[minus_1]*(-3*log(theta[update_theta]) + log(theta[minus_2]) - 2*log(theta[minus_1])) + 
                                                                    w_k[plus_1]*(-2*log(theta[plus_1]) + log(theta[plus_2]) - 3*log(theta[update_theta])))
            U_2 <- (16*w_k[update_theta] + 8*w_k[minus_1] + 8*w_k[plus_1]) * lambda
            U_3 <- temp4[-c(1,2,num_nonzero-1,num_nonzero)]
          }
          theta[non_zero_theta[1]] <- tryCatch(uniroot(g_1,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                               error = function(e) {return(theta[non_zero_theta[1]])})
          
          theta[non_zero_theta[2]] <- tryCatch(uniroot(g_2,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                               error = function(e) {return(theta[non_zero_theta[2]])})
          #parallelization here
          if (length(update_theta) > 0)
            for (jj in 1:length(update_theta)) {
              g <- function(x){U_1[jj]/x - U_2[jj]*log(x)/x - U_3[jj]}
              theta[update_theta[jj]] <- tryCatch(uniroot(g,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                  error = function(e) theta[update_theta[jj]])
            }
          theta[non_zero_theta[num_nonzero - 1]] <- tryCatch(uniroot(g_semiend,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                             error = function(e) 
                                                             {return(theta[non_zero_theta[num_nonzero - 1]])})
          theta[non_zero_theta[num_nonzero]] <- tryCatch(uniroot(g_end,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                         error = function(e) {return(theta[non_zero_theta[num_nonzero]])})
        }
        else if (1 == mode_reg_A) {
          #mode_reg_A == 1  
          g_1  <- function(x) {
            (net_stat$sum_m_k[non_zero_theta[ok_log10[1]]] - 2 * w_k[non_zero_theta[ok_log10[2]]] * lambda * (-2 * log(theta[non_zero_theta[ok_log10[2]]]) / 
                                                                                                                log(center_k[non_zero_theta[ok_log10[2]]])  + 
                                                                                                                log(theta[non_zero_theta[ok_log10[3]]]) / log(center_k[non_zero_theta[ok_log10[3]]]) - 
                                                                                                                3 * log(theta[non_zero_theta[ok_log10[1]]]) / log(center_k[non_zero_theta[ok_log10[1]]])  ))/x -   
              temp4[ok_log10[1]] - 8 * lambda * w_k[non_zero_theta[ok_log10[2]]] * log(x) / x / log(center_k[non_zero_theta[ok_log10[1]]])}
          
          if (length(ok_log10) > 3) {
            g_2 <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[ok_log10[2]]] - 2 * lambda * 
                 (2 * w_k[non_zero_theta[ok_log10[2]]] * (- 2 * log(theta[non_zero_theta[ok_log10[2]]]) / log(center_k[non_zero_theta[ok_log10[2]]]) - 
                                                            log(theta[non_zero_theta[ok_log10[3]]]) / log(center_k[non_zero_theta[ok_log10[3]]]) - 
                                                            log(theta[non_zero_theta[ok_log10[1]]]) / log(center_k[non_zero_theta[ok_log10[1]]])   ) + 
                    w_k[non_zero_theta[ok_log10[3]]]*(- 2 * log(theta[non_zero_theta[ok_log10[3]]]) / log(center_k[non_zero_theta[ok_log10[3]]])  +
                                                        log(theta[non_zero_theta[ok_log10[4]]]) / log(center_k[non_zero_theta[ok_log10[4]]]) -
                                                        3 * log(theta[non_zero_theta[ok_log10[2]]]) / log(center_k[non_zero_theta[ok_log10[2]]])  ))) / x - 
                temp4[ok_log10[2]] - (16*w_k[non_zero_theta[ok_log10[2]]] + 8*w_k[non_zero_theta[ok_log10[3]]]) * lambda * log(x) / x / log(center_k[non_zero_theta[ok_log10[2]]])}
          } else if (length(ok_log10) == 3) {
            g_2 <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[ok_log10[2]]] - 2 * lambda * 
                 (2 * w_k[non_zero_theta[ok_log10[2]]] * (- 2 * log(theta[non_zero_theta[ok_log10[2]]]) / log(center_k[non_zero_theta[ok_log10[2]]]) - 
                                                            log(theta[non_zero_theta[ok_log10[3]]]) / log(center_k[non_zero_theta[ok_log10[3]]]) - 
                                                            log(theta[non_zero_theta[ok_log10[1]]]) / log(center_k[non_zero_theta[ok_log10[1]]])   ) + 
                    w_k[non_zero_theta[ok_log10[3]]]*(- 2 * log(theta[non_zero_theta[ok_log10[3]]]) / log(center_k[non_zero_theta[ok_log10[3]]]) -
                                                        3 * log(theta[non_zero_theta[ok_log10[2]]]) / log(center_k[non_zero_theta[ok_log10[2]]])  ))) / x - 
                temp4[ok_log10[2]] - (16*w_k[non_zero_theta[ok_log10[2]]] + 8*w_k[non_zero_theta[ok_log10[3]]]) * lambda * log(x) / x / log(center_k[non_zero_theta[ok_log10[2]]])}
          }
          if (num_ok_log10 > 3) {
            g_semiend <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] - 2*lambda * 
                 (2*w_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] * ( -2 * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) - 
                                                                          log(theta[non_zero_theta[ok_log10[num_ok_log10]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10]]]) - 
                                                                          log(theta[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) ) +
                    w_k[non_zero_theta[ok_log10[num_ok_log10 -2]]] *( - 3 * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) + 
                                                                        log(theta[non_zero_theta[ok_log10[num_ok_log10 - 3]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 3]]]) - 
                                                                        2 * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) ) ) ) / x - 
                temp4[ok_log10[num_ok_log10 - 1]] - 
                (16*w_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] + 8 * w_k[non_zero_theta[ok_log10[num_ok_log10 -2]]]) * lambda * log(x) / x / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]])}
          } else if (num_ok_log10 == 3) {
            g_semiend <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] - 2*lambda * 
                 (2*w_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] * ( -2 * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) - 
                                                                          log(theta[non_zero_theta[ok_log10[num_ok_log10]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10]]]) - 
                                                                          log(theta[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) ) +
                    w_k[non_zero_theta[ok_log10[num_ok_log10 -2]]] *( - 3 * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) - 
                                                                        2 * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) ) ) ) / x - 
                temp4[ok_log10[num_ok_log10 - 1]] - 
                (16*w_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] + 8 * w_k[non_zero_theta[ok_log10[num_ok_log10 -2]]]) * lambda * log(x) / x / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]])}   
          }
          
          g_end     <-  function(x) {
            (net_stat$sum_m_k[non_zero_theta[ok_log10[num_ok_log10]]] - 2 * lambda * w_k[non_zero_theta[ok_log10[num_ok_log10 -1]]] * (-3 * 
                                                                                                                                         log(theta[non_zero_theta[ok_log10[num_ok_log10]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10]]]) + 
                                                                                                                                         log(theta[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) - 
                                                                                                                                         2 * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) / log(center_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]])   )) / x  -   
              temp4[ok_log10[num_ok_log10]] - 8 * w_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] * lambda * log(x)/x / log(center_k[non_zero_theta[ok_log10[num_ok_log10]]]) }  
          
          
          if (length(update_theta) > 0) {
            U_1 <- net_stat$sum_m_k[update_theta] - 2 * lambda * (2 * w_k[update_theta] * (-2 * log(theta[update_theta]) / log(center_k[update_theta]) - 
                                                                                             log(theta[plus_1]) /  log(center_k[plus_1]) - log(theta[minus_1]) / log(center_k[minus_1]) ) + 
                                                                    w_k[minus_1]*(-3 * log(theta[update_theta]) / log(center_k[update_theta]) + 
                                                                                    log(theta[minus_2]) / log(center_k[minus_2]) - 2 * log(theta[minus_1]) / log(center_k[minus_1])  ) + 
                                                                    w_k[plus_1]*(-2 * log(theta[plus_1]) / log(center_k[plus_1]) + log(theta[plus_2]) / log(center_k[plus_2]) - 
                                                                                   3 * log(theta[update_theta]) / log(center_k[update_theta])))
            
            U_2 <- (16*w_k[update_theta] + 8*w_k[minus_1] + 8*w_k[plus_1]) * lambda
            
            U_3 <- temp4[ok_log10[-c(1,2,num_ok_log10 - 1,num_ok_log10)]]
            
          }
          
          
          theta[non_zero_theta[ok_log10[1]]] <- tryCatch(uniroot(g_1,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                         error = function(e) {return(theta[non_zero_theta[ok_log10[1]]])})
          
          theta[non_zero_theta[ok_log10[2]]] <- tryCatch(uniroot(g_2,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                         error = function(e) {return(theta[non_zero_theta[ok_log10[2]]])})
          #parallelization here
          if (length(update_theta) > 0)
            for (jj in 1:length(update_theta)) {
              g <- function(x){U_1[jj]/x - U_2[jj]*log(x)/x / log(center_k[update_theta[jj]]) - U_3[jj]}
              theta[update_theta[jj]] <- tryCatch(uniroot(g,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                  error = function(e) theta[update_theta[jj]])
            }
          theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]] <- tryCatch(uniroot(g_semiend,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                                        error = function(e) 
                                                                        {return(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]])})
          theta[non_zero_theta[ok_log10[num_ok_log10]]] <- tryCatch(uniroot(g_end,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                                    error = function(e) {return(theta[non_zero_theta[ok_log10[num_ok_log10]]])})
          
          if (length(not_ok_log10) > 0) 
            theta[non_zero_theta[not_ok_log10]] <- net_stat$sum_m_k[non_zero_theta[not_ok_log10]]/temp4[not_ok_log10]  
        } else if (2 == mode_reg_A) {
          #mode_reg_A == 2  
          
          g_1_new <- function(x) {
            (net_stat$sum_m_k[non_zero_theta[ok_log10[1]]] - 2 * u[2,1] * w_k[non_zero_theta[ok_log10[2]]] * lambda * 
               (- u[2,2] * log(theta[non_zero_theta[ok_log10[2]]])   + 
                  u[2,3] * log(theta[non_zero_theta[ok_log10[3]]])  - 
                  3 * u[2,1] * log(theta[non_zero_theta[ok_log10[1]]])  ) + 
               2 * lambda * w_k[non_zero_theta[not_ok_log10]] * 
               (log(theta[non_zero_theta[ok_log10[1]]])  + log(theta[non_zero_theta[not_ok_log10]]) )) / x -   
              temp4[ok_log10[1]] - (8 * lambda * w_k[non_zero_theta[ok_log10[2]]] * u[2,1]^2 + 
                                      4 * lambda * w_k[non_zero_theta[not_ok_log10]])* log(x) / x }
          if (length(ok_log10) > 3) {
            g_2 <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[ok_log10[2]]] - lambda * 
                 (-2 * u[2 , 2] * w_k[non_zero_theta[ok_log10[2]]] * ( u[2 , 2] * log(theta[non_zero_theta[ok_log10[2]]]) + 
                                                                         u[2 , 3] * log(theta[non_zero_theta[ok_log10[3]]]) +  
                                                                         u[2 , 1] * log(theta[non_zero_theta[ok_log10[1]]]) ) +
                    
                    2 * u[3,2] * w_k[non_zero_theta[ok_log10[3]]]*(- u[3 , 3] * log(theta[non_zero_theta[ok_log10[3]]])  +
                                                                     u[3 , 4] * log(theta[non_zero_theta[ok_log10[4]]]) - 
                                                                     3 * u[3 , 2] * log(theta[non_zero_theta[ok_log10[2]]])  ) )) / x - 
                
                temp4[ok_log10[2]] + (-4 * u[2,2]^2 * w_k[non_zero_theta[ok_log10[2]]] - 
                                        8 * u[3,2]^2 * w_k[non_zero_theta[ok_log10[3]]]) * lambda * log(x) / x     }
          } else if (length(ok_log10) == 3) {
            g_2 <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[ok_log10[2]]] - lambda * 
                 (-2 * u[2 , 2] * w_k[non_zero_theta[ok_log10[2]]] * ( u[2 , 2] * log(theta[non_zero_theta[ok_log10[2]]]) + 
                                                                         u[2 , 3] * log(theta[non_zero_theta[ok_log10[3]]]) +  
                                                                         u[2 , 1] * log(theta[non_zero_theta[ok_log10[1]]]) ) +
                    
                    2 * u[3,2] * w_k[non_zero_theta[ok_log10[3]]]*(- u[3 , 3] * log(theta[non_zero_theta[ok_log10[3]]]) - 
                                                                     3 * u[3 , 2] * log(theta[non_zero_theta[ok_log10[2]]])  ) )) / x - 
                
                temp4[ok_log10[2]] + (-4 * u[2,2]^2 * w_k[non_zero_theta[ok_log10[2]]] - 
                                        8 * u[3,2]^2 * w_k[non_zero_theta[ok_log10[3]]]) * lambda * log(x) / x     }  
            
          }
          if (num_ok_log10 > 3) {
            g_semiend <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] - lambda * 
                 (-2 * u[num_ok_log10 - 1, num_ok_log10 - 1] * w_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] * 
                    ( u[num_ok_log10 - 1, num_ok_log10 - 1] * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) + 
                        u[num_ok_log10 - 1, num_ok_log10] * log(theta[non_zero_theta[ok_log10[num_ok_log10]]])  + 
                        u[num_ok_log10 - 1, num_ok_log10 -  2] * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) )  +
                    2 * u[num_ok_log10 - 2, num_ok_log10 - 1] * w_k[non_zero_theta[ok_log10[num_ok_log10 -2]]] * 
                    ( - 3 * u[num_ok_log10 - 2, num_ok_log10 - 1]* log(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) + 
                        u[num_ok_log10 - 2, num_ok_log10 - 3] * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 3]]])  - 
                        u[num_ok_log10 - 2, num_ok_log10 - 2] * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 2]]])  ) ) ) / x - 
                temp4[ok_log10[num_ok_log10 - 1]] - 
                (4 * u[num_ok_log10 - 1, num_ok_log10 - 1]^2 * w_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] + 8 * u[num_ok_log10 - 2, num_ok_log10 - 1] ^ 2 * 
                   w_k[non_zero_theta[ok_log10[num_ok_log10 -2]]]) * lambda * log(x) / x }
          } else if (num_ok_log10 == 3) {
            g_semiend <- function(x) {
              (net_stat$sum_m_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] - lambda * 
                 (-2 * u[num_ok_log10 - 1, num_ok_log10 - 1] * w_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] * 
                    (u[num_ok_log10 - 1, num_ok_log10 - 1] * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) + 
                       u[num_ok_log10 - 1, num_ok_log10] * log(theta[non_zero_theta[ok_log10[num_ok_log10]]])  + 
                       u[num_ok_log10 - 1, num_ok_log10 -  2] * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 2]]]) )  +
                    2 * u[num_ok_log10 - 2, num_ok_log10 - 1] * w_k[non_zero_theta[ok_log10[num_ok_log10 -2]]] * 
                    ( - 3 * u[num_ok_log10 - 2, num_ok_log10 - 1]* log(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]]) - 
                        u[num_ok_log10 - 2, num_ok_log10 - 2] * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 2]]])  ) ) ) / x - 
                temp4[ok_log10[num_ok_log10 - 1]] - 
                (4 * u[num_ok_log10 - 1, num_ok_log10 - 1]^2 * w_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] + 8 * u[num_ok_log10 - 2, num_ok_log10 - 1] ^ 2 * 
                   w_k[non_zero_theta[ok_log10[num_ok_log10 -2]]]) * lambda * log(x) / x }  
            
          }
          
          g_end     <-  function(x) {
            (net_stat$sum_m_k[non_zero_theta[ok_log10[num_ok_log10]]] - 
               2 * u [num_ok_log10 - 1, num_ok_log10] * lambda * w_k[non_zero_theta[ok_log10[num_ok_log10 -1]]] * 
               (-3 * u [num_ok_log10 - 1, num_ok_log10] * log(theta[non_zero_theta[ok_log10[num_ok_log10]]])  + 
                  u [num_ok_log10 - 1, num_ok_log10 - 2 ] * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 2]]])  - 
                  u [num_ok_log10 - 1, num_ok_log10 - 1] * log(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]])   )) / x  -   
              temp4[ok_log10[num_ok_log10]] - 
              8 * u [num_ok_log10 - 1, num_ok_log10]^2 * w_k[non_zero_theta[ok_log10[num_ok_log10 - 1]]] * lambda * log(x) / x }  
          
          if (length(update_theta) > 0) {
            U_1 <- net_stat$sum_m_k[update_theta] - lambda * (- 2 * extract_u(update_k, update_k) * w_k[update_theta] * (
              extract_u(update_k, update_k) * log(theta[update_theta]) + 
                extract_u(update_k, update_k + 1) * log(theta[plus_1]) +
                extract_u(update_k, update_k - 1) * log(theta[minus_1]) ) +
                
                2 * extract_u(update_k - 1, update_k) * w_k[minus_1] * (
                  -3 * extract_u(update_k - 1, update_k) * log(theta[update_theta]) + 
                    extract_u(update_k - 1, update_k - 2) * log(theta[minus_2]) - 
                    extract_u(update_k - 1, update_k - 1) * log(theta[minus_1]) ) + 
                
                2 * extract_u(update_k + 1, update_k) * w_k[plus_1]* (
                  - extract_u(update_k + 1, update_k + 1) * log(theta[plus_1]) + 
                    extract_u(update_k + 1, update_k + 2) * log(theta[plus_2])  - 
                    3 * extract_u(update_k + 1, update_k)* log(theta[update_theta]) )
            )
            
            U_2 <- (4 * extract_u(update_k, update_k) ^ 2 * w_k[update_theta] + 8 * extract_u(update_k - 1, update_k) ^ 2 * w_k[minus_1] + 
                      8 * extract_u(update_k + 1, update_k) ^ 2 * w_k[plus_1]) * lambda
            U_3 <- temp4[ok_log10[-c(1 , 2 , num_ok_log10 - 1 , num_ok_log10)]]
            
          }
          
          theta[non_zero_theta[ok_log10[1]]] <- tryCatch(uniroot(g_1_new,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                         error = function(e) {return(theta[non_zero_theta[ok_log10[1]]])})
          
          theta[non_zero_theta[ok_log10[2]]] <- tryCatch(uniroot(g_2,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                         error = function(e) {return(theta[non_zero_theta[ok_log10[2]]])})
          #parallelization here
          if (length(update_theta) > 0)
            for (jj in 1:length(update_theta)) {
              #jj <- length(update_theta)
              g <- function(x){U_1[jj]/x - U_2[jj]*log(x)/x - U_3[jj]}
              theta[update_theta[jj]] <- tryCatch(uniroot(g,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                  error = function(e) theta[update_theta[jj]])
            }
          theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]] <- tryCatch(uniroot(g_semiend,interval = c(0.0000001 , 1000),tol = .Machine$double.eps)$root,
                                                                        error = function(e) 
                                                                        {return(theta[non_zero_theta[ok_log10[num_ok_log10 - 1]]])})
          theta[non_zero_theta[ok_log10[num_ok_log10]]] <- tryCatch(uniroot(g_end,interval = c(0.0000001 , 1000),tol = .Machine$double.eps)$root,
                                                                    error = function(e) {return(theta[non_zero_theta[ok_log10[num_ok_log10]]])})
          
          if (length(not_ok_log10) > 0) {
           
            if (lambda <= 0) {
              theta[non_zero_theta[not_ok_log10]] <- net_stat$sum_m_k[non_zero_theta[not_ok_log10]]/temp4[not_ok_log10] 
            }else {
              g_0  <- function(x) {
                (net_stat$sum_m_k[non_zero_theta[not_ok_log10]] - 2 * w_k[non_zero_theta[not_ok_log10]] * lambda * 
                   (- log(theta[non_zero_theta[not_ok_log10]]) - 
                      log(theta[non_zero_theta[ok_log10[1]]]) )) / x -   
                  temp4[not_ok_log10] - 4 * lambda * w_k[non_zero_theta[not_ok_log10]] * log(x) / x}
              theta[non_zero_theta[not_ok_log10]] <- tryCatch(uniroot(g_0,interval = c(0.0000001 , 1000),tol = .Machine$double.eps)$root,
                                                              error = function(e) {return(theta[non_zero_theta[not_ok_log10]])}) 
            }
            
          }
          
        }
       
        
        non_zero_theta_now <- theta > 10^-30
        
        if (length(non_zero_theta_now) > 0)
          for (lll in 1:length(non_zero_theta_now))
            if (non_zero_theta_start[lll] != non_zero_theta_now[lll]) {
              
              diverge_zero_theta <- TRUE
              break; 
            }
        if (TRUE == diverge_zero_theta)
          break;   
      }
    } else {
      # log-linear PA
     
      #if (i == 1 || only_f == FALSE)
        .normalized_constant_alpha(normalized_const, alpha,
                                   PA_offset,net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
     
      time_non_zero     <- which(normalized_const != 0)
      #non_zero_f        <- which(f > 0) 
     
      
      upper_f_term <- net_stat$z_j[non_zero_f] * log(f[non_zero_f])
  
      
      log10_likelihood    <- c(log10_likelihood, sum(upper_f_term) +
                                 alpha * sum(net_stat$sum_m_k * log(theta)) -
                                 sum(net_stat$m_t[time_non_zero] * log(normalized_const[time_non_zero])) +
                                 (sum((shape / weight_f[non_zero_f] - 1) * log(f[non_zero_f])) - 
                                    sum(rate / weight_f[non_zero_f]  * (f[non_zero_f]))) +
                                 sum(net_stat$offset_m_tk)*log(offset) + (shape - 1) * log(offset) - rate * offset)
      
      
      log10_likelihood[length(log10_likelihood)] <- log10_likelihood[length(log10_likelihood)] + net_stat$sum_m_k[1] * 
                                                    log(PA_offset);
     
      names(log10_likelihood) <- NULL
      
      if ((TRUE == debug) && (length(log10_likelihood) > 0)){
        print(log10_likelihood[length(log10_likelihood)])
      }
      break_flag <- FALSE
      
      if (length(log10_likelihood) > 1) {
        if (abs(log10_likelihood[length(log10_likelihood)] - log10_likelihood[length(log10_likelihood) - 1]) /
            abs(log10_likelihood[length(log10_likelihood)] + 1) > 10^-15)  
          if (log10_likelihood[length(log10_likelihood)] < log10_likelihood[length(log10_likelihood) - 1]) {
           
            break_flag <- TRUE;
            
          }
      }
     
      
    
      if (TRUE == auto_stop)
        if (length(log10_likelihood) > 1)
            tryCatch({if (abs(log10_likelihood[length(log10_likelihood)] - log10_likelihood[length(log10_likelihood) - 1]) / 
                         (abs(log10_likelihood[length(log10_likelihood) - 1]) + 1) < stop_cond)
              break_flag <- TRUE;},error = function(e) { 
              break_flag <- TRUE;})
      
      
      
      if (break_flag)
        break;   
    
      # Remember here theta is always the degree sequence,i.e. we need to calculate the power alpha of this deg. seq.    
      
      ##################### Update f ######################
      #.normalized_constant_alpha(normalized_const, alpha,PA_offset, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
      if (FALSE == only_PA)  {
           
          .update_f_alpha_new(f,non_zero_f,alpha,PA_offset,net_stat$node_degree,theta,
                              net_stat$z_j,normalized_const,net_stat$m_t,shape,rate,
                              weight_f)
          .normalized_constant_alpha(normalized_const, alpha,PA_offset, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
        # update offset
        # if (net_stat$deg_thresh > 0)
        # offset <- .update_offset_alpha(alpha,net_stat$offset_tk, net_stat$offset_m_tk, theta, normalized_const,net_stat$m_t, shape,rate)
      }
      non_zero_f_now <- f > 10^-30
      
      if (length(non_zero_f_now) > 0)
        for (lll in 1:length(non_zero_f_now))
          if (non_zero_f_start[lll] != non_zero_f_now[lll]) {
            diverge_zero <- TRUE
            break; 
          }
      if (TRUE == diverge_zero)
        break;   
      
      #if (sum(f) != 0)            
      #     if(normalized_f == TRUE) 
      #         f <- f / sum(f) * length(f)   
      #####################  Update alpha #######################  
      
     
      if (FALSE == only_f) 
        alpha     <- .update_alpha_fast(non_zero_theta,
                                        normalized_const,f, PA_offset,
                                        theta , net_stat$node_degree,net_stat$m_t,net_stat$sum_m_k,net_stat$offset_tk,
                                        offset, alpha) 
     #.normalized_constant_alpha(normalized_const, alpha,PA_offset, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
      #PA_offset <- .update_PA_offset(normalized_const,f,net_stat$node_degree,net_stat$m_t,net_stat$sum_m_k,
      #                                 net_stat$offset_tk);
      #.normalized_constant_alpha(normalized_const, alpha, PA_offset, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
      alpha_series <- c(alpha_series,alpha)
      
      theta_temp         <- theta^alpha
      non_zero_theta_now <- theta_temp > 10^-30
      
      if (length(non_zero_theta_now) > 0)
        for (lll in 1:length(non_zero_theta_now))
          if (non_zero_theta_start[lll] != non_zero_theta_now[lll]) {
            diverge_zero_theta <- TRUE
            break; 
          }
      if (TRUE == diverge_zero_theta)
        break; 
      
    }
  }
  
  if (is.null(true_A)) {
      if (mode_f == "Log_linear")   
          theta    <- theta^alpha  
      if (only_f == TRUE)
          theta[1] <- PA_offset  
  }
  
  ############################################   ############################################   
  ############################################   ############################################   
  ########## ########################### End of Iteration ############################################   
  ############################################   ############################################   
  ############################################   ############################################   
  
   
 
  if (normalized_f == TRUE) {
    sum_f  <- sum(f)  
    f      <- (length(f) + 1) * f / (sum_f + offset)
    offset <- offset * (length(f) + 1) / (sum_f + offset)   
    if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))){
      .normalized_constant(normalized_const,net_stat$node_degree,theta,f,net_stat$offset_tk,offset) 
    }
    else normalized_const <- as.vector(net_stat$n_tk[,non_zero_theta]%*%theta[non_zero_theta])
  }
  
  
  time_non_zero       <- which(normalized_const > 0)
  non_zero_f_temp     <- which(f > 0)
  non_zero_theta_temp <- which(theta > 0) 
  
    
  
  
  ##### Variance of f ####################
  cov_f   <- rep(0,length(f))
  if (FALSE == only_PA)
    .cal_var_f_new(cov_f,non_zero_f,net_stat$node_degree,theta,f,net_stat$z_j,
               normalized_const,net_stat$m_t,shape,weight_f)
  

  
  ### Variance of alpha if mode_f == "Log_linear" ######
  if (mode_f == "Log_linear") {
    .normalized_constant_alpha(normalized_const, alpha,PA_offset, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
    cov_alpha <- .var_alpha(alpha, non_zero_theta,
                            normalized_const,f, PA_offset,
                            theta,net_stat$node_degree,net_stat$m_t,net_stat$sum_m_k,net_stat$offset_tk,offset) 
  }
  

  # mode_reg_A = 0:
  hessian_of_regularization <- function(theta){
    n      <- length(theta)
    result <- vector()
    if (n - 2 >= 3){
      result <- c(result,2*w_k[non_zero_theta[2]]*((1-log(theta[1]))/theta[1]^2 + (2*log(theta[2]) - log(theta[3])) / theta[1]^2 ) )
      
      result <- c(result,2*w_k[non_zero_theta[3]]*((1-log(theta[2]))/theta[2]^2 + (2*log(theta[3]) - log(theta[4])) / theta[2]^2 ) +
                    2*w_k[non_zero_theta[2]]*((2-2*log(theta[2]))/theta[2]^2 + (log(theta[3]) + log(theta[1])) / theta[2]^2 ))
      for (ii in 3:(n-2))
        result <- c(result,2*w_k[non_zero_theta[ii]]*((2-2*log(theta[ii]))/theta[ii]^2 + (log(theta[ii+1]) + log(theta[ii-1])) / theta[ii]^2 ) + 
                      2*w_k[non_zero_theta[ii + 1]]*((1-log(theta[ii]))/theta[ii]^2 + (2*log(theta[ii+1]) - log(theta[ii+2])) / theta[ii]^2 ) +
                      2*w_k[non_zero_theta[ii-1]]*((1-log(theta[ii]))/theta[ii]^2 + (2*log(theta[ii-1]) - log(theta[ii-2])) / theta[ii]^2 ))
      ii <- n - 1
      result <- c(result, 2*w_k[non_zero_theta[ii]]*((2-2*log(theta[ii]))/theta[ii]^2 + (log(theta[ii+1]) + log(theta[ii-1])) / theta[ii]^2 ) +
                    2*w_k[non_zero_theta[ii-1]]*((1-log(theta[ii]))/theta[ii]^2 + (2*log(theta[ii-1]) - log(theta[ii-2])) / theta[ii]^2 ))
      ii <- n
      result <- c(result, 2*w_k[non_zero_theta[ii-1]]*((1-log(theta[ii]))/theta[ii]^2 + (2*log(theta[ii-1]) - log(theta[ii-2])) / theta[ii]^2 )) 
    }
    return(w_k[non_zero_theta]*result)
  }
  
  # mode_reg_A = 1:
  hessian_of_regularization_log10 <- function(theta){
    n      <- length(theta)
    result <- vector()
    
    result <- c(result,2*w_k[non_zero_theta[ok_log10[2]]]*((1 - log(theta[1]) / log(center_k[non_zero_theta[ok_log10[1]]]) ) /theta[1]^2 + 
                                                             (2*log(theta[2]) / log(center_k[non_zero_theta[ok_log10[2]]]) - 
                                                                log(theta[3])) / log(center_k[non_zero_theta[ok_log10[3]]]) / theta[1]^2 ) )
    
    result <- c(result,2*w_k[non_zero_theta[ok_log10[3]]]*((1 - log(theta[2]) / log(center_k[non_zero_theta[ok_log10[2]]]) )/theta[2]^2 + 
                                                             (2*log(theta[3]) / log(center_k[non_zero_theta[ok_log10[3]]]) - log(theta[4]) / log(center_k[non_zero_theta[ok_log10[4]]])) / theta[2]^2 ) +
                  2*w_k[non_zero_theta[2]]*((2 - 2 * log(theta[2]) / log(center_k[non_zero_theta[ok_log10[2]]]) )/theta[2]^2 + ( log(theta[3]) / log(center_k[non_zero_theta[ok_log10[3]]]) + 
                                                                                                                                   log(theta[1]) / log(center_k[non_zero_theta[ok_log10[1]]])  ) / theta[2]^2 ))
    if (n -2 >= 3) { 
      for (ii in 3:(n-2))
        result <- c(result,2*w_k[non_zero_theta[ok_log10[ii]]]*((2 - 2 * log(theta[ii]) / log(center_k[non_zero_theta[ok_log10[ii]]]) )/theta[ii]^2 + 
                                                                  (log(theta[ii + 1]) / log(center_k[non_zero_theta[ok_log10[ii + 1]]]) + 
                                                                     log(theta[ii - 1]) / log(center_k[non_zero_theta[ok_log10[ii - 1]]]) ) / theta[ii]^2 ) + 
                      2*w_k[non_zero_theta[ok_log10[ii + 1]]]*((1 - log(theta[ii]) / log(center_k[non_zero_theta[ok_log10[ii]]]) )/theta[ii]^2 + 
                                                                 (2 * log(theta[ii + 1] / log(center_k[non_zero_theta[ok_log10[ii + 1]]])) - 
                                                                    log(theta[ii + 2]) / log(center_k[non_zero_theta[ok_log10[ii + 2]]]) ) / theta[ii]^2 ) +
                      2*w_k[non_zero_theta[ok_log10[ii-1]]]*((1 - log(theta[ii]) / log(center_k[non_zero_theta[ok_log10[ii]]]) )/theta[ii]^2 + 
                                                               (2 * log(theta[ii - 1]) / log(center_k[non_zero_theta[ok_log10[ii - 1]]]) - 
                                                                  log(theta[ii - 2]) / log(center_k[non_zero_theta[ok_log10[ii - 2]]]) ) / theta[ii]^2 ))
    }
    ii <- n - 1
    result <- c(result, 2 * w_k[non_zero_theta[ok_log10[ii]]]*((2 - 2 * log(theta[ii]) / log(center_k[non_zero_theta[ok_log10[ii]]]) )/theta[ii]^2 + 
                                                                 (log(theta[ii + 1]) / log(center_k[non_zero_theta[ok_log10[ii + 1]]]) + 
                                                                    log(theta[ii - 1]) / log(center_k[non_zero_theta[ok_log10[ii - 1]]])  ) / theta[ii]^2 ) +
                  2*w_k[non_zero_theta[ok_log10[ii-1]]]*((1 - log(theta[ii]) / log(center_k[non_zero_theta[ok_log10[ii]]])) / theta[ii]^2 + 
                                                           (2 * log(theta[ii - 1]) / log(center_k[non_zero_theta[ok_log10[ii - 1]]]) - 
                                                              log(theta[ii - 2]) / log(center_k[non_zero_theta[ok_log10[ii - 2]]]) ) / theta[ii]^2 ))
    ii <- n
    result <- c(result, 2*w_k[non_zero_theta[ok_log10[ii-1]]]*((1 - log(theta[ii]) / log(center_k[non_zero_theta[ok_log10[ii]]]) )/theta[ii]^2 + 
                                                                 (2 * log(theta[ii - 1]) / log(center_k[non_zero_theta[ok_log10[ii - 1]]]) - 
                                                                    log(theta[ii - 2]) / log(center_k[non_zero_theta[ok_log10[ii - 2]]])  ) / theta[ii]^2 )) 
    if (length(result) == 0)
      result <- rep(1,length(w_k[non_zero_theta[ok_log10]]))
    return(w_k[non_zero_theta[ok_log10]] * result)
  }
  
  # mode_reg_A = 2:
  hessian_of_regularization_mode_2 <- function(theta){
    n      <- length(theta)
    result <- vector()
    if ((n >= 2) && (num_ok_log10 >= 3)) {
      if (length(not_ok_log10) > 0) {  
        result <- c(result, 2 * w_k[non_zero_theta[not_ok_log10]] * (1 - log(theta[non_zero_theta[1]]) + log(theta[non_zero_theta[ok_log10[1]]])) / theta[non_zero_theta[1]]^2) 
        
        result <- c(result, - (-2 *  w_k[non_zero_theta[not_ok_log10]] * log(theta[non_zero_theta[not_ok_log10]]) - 2 *  w_k[non_zero_theta[not_ok_log10]]  + 
                                 2 *  w_k[non_zero_theta[not_ok_log10]] * log(theta[non_zero_theta[ok_log10[1]]])  +  
                                 2 *  w_k[non_zero_theta[ok_log10[2]]] * u[2,1] * (u[2,3] * log(theta[non_zero_theta[ok_log10[3]]]) -  u[2,2] * log(theta[non_zero_theta[ok_log10[2]]])) - 
                                 2 *  w_k[non_zero_theta[ok_log10[2]]] * u[2,1]^2 + 2 *  w_k[non_zero_theta[ok_log10[2]]] * u[2,1]^2 * log(theta[non_zero_theta[ok_log10[1]]]) ) / theta[non_zero_theta[ok_log10[1]]]^2);
      } else {
        result <- c(result, 0) 
        
        result <- c(result, - (  2 *  w_k[non_zero_theta[not_ok_log10]] * log(theta[non_zero_theta[ok_log10[1]]])  +  
                                   2 *  w_k[non_zero_theta[ok_log10[2]]] * u[2,1] * (u[2,3] * log(theta[non_zero_theta[ok_log10[3]]]) -  u[2,2] * log(theta[non_zero_theta[ok_log10[2]]])) - 
                                   2 *  w_k[non_zero_theta[ok_log10[2]]] * u[2,1]^2 + 2 *  w_k[non_zero_theta[ok_log10[2]]] * u[2,1] ^ 2 * log(theta[non_zero_theta[ok_log10[1]]]) ) / theta[non_zero_theta[ok_log10[1]]] ^ 2); 
      }
      
      # result <- c(result, - (-2 *  w_k[non_zero_theta[ok_log10[2]]] * (u[2 , 3] * log(theta[non_zero_theta[ok_log10[3]]]) +  u[2,1] * log(theta[non_zero_theta[ok_log10[1]]]))  * u[2,2] - 
      #                         2 *  w_k[non_zero_theta[ok_log10[2]]] * u[2 , 2] ^ 2 + 
      #                         2 *  w_k[non_zero_theta[ok_log10[2]]] * u[2 , 2] ^ 2 *  log(theta[non_zero_theta[ok_log10[2]]]) + 
      #                         2 *  w_k[non_zero_theta[ok_log10[3]]] * (u[3 , 4] * log(theta[non_zero_theta[ok_log10[4]]]) -  u[3,3] * log(theta[non_zero_theta[ok_log10[3]]]))  * u[3,2] - 
      #                         2 *  w_k[non_zero_theta[ok_log10[3]]] * u[3 , 2]^2 +  2 *  w_k[non_zero_theta[ok_log10[3]]] * u[3,2]^2 * log(theta[non_zero_theta[ok_log10[2]]]) ) / theta[non_zero_theta[ok_log10[2]]] ^ 2);
      result <- c(result, - (-2 *  w_k[non_zero_theta[ok_log10[2]]] * (u[2 , 3] * log(theta[non_zero_theta[ok_log10[3]]]) +  u[2,1] * log(theta[non_zero_theta[ok_log10[1]]]))  * u[2,2] - 
                               2 *  w_k[non_zero_theta[ok_log10[2]]] * u[2 , 2] ^ 2 + 
                               2 *  w_k[non_zero_theta[ok_log10[2]]] * u[2 , 2] ^ 2 *  log(theta[non_zero_theta[ok_log10[2]]]) + 
                               2 *  w_k[non_zero_theta[ok_log10[3]]] * ( -  u[3,3] * log(theta[non_zero_theta[ok_log10[3]]]))  * u[3,2] - 
                               2 *  w_k[non_zero_theta[ok_log10[3]]] * u[3 , 2]^2 +  2 *  w_k[non_zero_theta[ok_log10[3]]] * u[3,2]^2 * log(theta[non_zero_theta[ok_log10[2]]]) ) / theta[non_zero_theta[ok_log10[2]]] ^ 2); 
      if (dim(u)[1] >= 4)
        result[length(result)] <- result[length(result)] - 2 *  w_k[non_zero_theta[ok_log10[3]]] * (u[3 , 4] * log(theta[non_zero_theta[ok_log10[4]]])) * u[3,2] / theta[non_zero_theta[ok_log10[2]]] ^ 2 
      m <- length(ok_log10)
      if (m - 2 >= 3)
        for (ii in 3:(m-2))
          result <- c(result, - ( 2 *  w_k[non_zero_theta[ok_log10[ii - 1]]] * ( - u[ii - 1,ii] + log(theta[non_zero_theta[ok_log10[ii]]]) + u[ii - 1,ii - 2] * log(theta[non_zero_theta[ok_log10[ii - 2]]]) -
                                                                                   u[ii - 1, ii - 1] * log(theta[non_zero_theta[ok_log10[ii - 1]]]))  * u[ii - 1,ii]  + 
                                    2 *  w_k[non_zero_theta[ok_log10[ii + 1]]] * ( u[ii + 1,ii + 2] *  log(theta[non_zero_theta[ok_log10[ii + 2]]]) -  u[ii + 1,ii] + log(theta[non_zero_theta[ok_log10[ii]]]) -
                                                                                     u[ii + 1,ii + 1] * log(theta[non_zero_theta[ok_log10[ii + 1]]]))  * u[ii + 1,ii] -    
                                    2 *  w_k[non_zero_theta[ok_log10[ii]]] * ( u[ii,ii + 1] *  log(theta[non_zero_theta[ok_log10[ii + 1]]]) +  u[ii,ii] - log(theta[non_zero_theta[ok_log10[ii]]]) +
                                                                                 u[ii,ii - 1] * log(theta[non_zero_theta[ok_log10[ii - 1]]]))  * u[ii,ii]  ) / theta[non_zero_theta[ok_log10[ii]]]^2)
      
      ii <- m - 1
      result <- c(result, - ( - 2 *  w_k[non_zero_theta[ok_log10[ii]]] * (u[ii,ii + 1] * log(theta[non_zero_theta[ok_log10[ii + 1]]]) + u[ii,ii - 1] * log(theta[non_zero_theta[ok_log10[ii - 1]]]) +
                                                                            u[ii, ii] - log(theta[non_zero_theta[ok_log10[ii]]]))  * u[ii,ii]  + 
                                2 *  w_k[non_zero_theta[ok_log10[ii - 1]]] * ( - u[ii - 1, ii]  +  log(theta[non_zero_theta[ok_log10[ii]]]) + u[ii - 1,ii - 2] *  log(theta[non_zero_theta[ok_log10[ii - 2]]]) -  
                                                                                 u[ii - 1,ii - 1] * log(theta[non_zero_theta[ok_log10[ii - 1]]]) ) * u[ii - 1,ii]  ) / theta[non_zero_theta[ok_log10[ii]]]^2) 
      ii <- m
      result <- c(result, - (2 *  w_k[non_zero_theta[ok_log10[ii - 1]]] * (- u[ii - 1, ii] +  log(theta[non_zero_theta[ok_log10[ii]]]) + u[ii - 1,ii - 2] * log(theta[non_zero_theta[ok_log10[ii - 2]]]) - 
                                                                             u[ii - 1,ii - 1] * log(theta[non_zero_theta[ok_log10[ii - 1]]]) ) * u[ii - 1,ii] / theta[non_zero_theta[ok_log10[ii]]]^2 )) 
    }
    return(result);
    
  }
  
  #interpolation for PA function
  #but only if (FALSE == only_f)
  if ((FALSE == only_f) && (mode_f[1] != "Log_linear")) {  
    if (TRUE == interpolate) {
      theta_nonzero <- which(net_stat$sum_m_k != 0)
      if ((only_f == TRUE) && ((mode_f[1] = "Linear_PA") || (mode_f[1] = "Log_linear"))) {
        if (0 != PA_offset)    
          theta_nonzero <- c(1,theta_nonzero)  
      }
      if (length(theta_nonzero) > 0)
        if (theta_nonzero[1] > 1) {
          theta[1:(theta_nonzero[1] - 1)]   <- theta[theta_nonzero[1]]  
          #cov_bin[1:(theta_nonzero[1] - 1)] <- cov_bin[theta_nonzero[1]] 
        }
      if (length(theta_nonzero) > 1) {
        for (i in 1:(length(theta_nonzero) - 1))
          if (theta_nonzero[i + 1] > theta_nonzero[i] + 1) 
            if (center_k[theta_nonzero[i]] > 0 && center_k[theta_nonzero[i + 1]] > 0 &&
                theta[theta_nonzero[i]] > 0 && theta[theta_nonzero[i + 1]] > 0)  {
              regress <- lm(c(log(theta[theta_nonzero[i]]) , log(theta[theta_nonzero[i + 1]]))~ 
                              c(log(center_k[theta_nonzero[i]]) , log(center_k[theta_nonzero[i + 1]]))) 
              for (j in (theta_nonzero[i] + 1) : (theta_nonzero[i + 1] - 1))
                theta[j] <- exp(log(center_k[j]) * regress$coefficients[2] + regress$coefficients[1])           
            }  
      }
      if (length(theta_nonzero) > 0)
        if (theta_nonzero[length(theta_nonzero)] < length(theta))
          theta[(theta_nonzero[length(theta_nonzero)] + 1):length(theta)] <- 
            theta[theta_nonzero[length(theta_nonzero)]]  
    }
  }
  
  beg     <- which(center_k >= net_stat$deg_thresh & theta != 0)[1]
  if (length(beg) == 0) beg <- which(theta != 0)[1]    
  theta <- theta/theta[beg]
  
  
  if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))){
    .normalized_constant(normalized_const,net_stat$node_degree,theta,f,net_stat$offset_tk,offset) 
  }
  else normalized_const <- as.vector(net_stat$n_tk[,non_zero_theta]%*%theta[non_zero_theta])
  ############# Calculate the variance of theta #####################
  cov_bin1       <- rep(0,length(theta)) 
  #non_zero       <- which(theta != 0)
  non_zero <- which(net_stat$sum_m_k > 0)
  non_zero_theta <- non_zero
  time_non_zero     <- which(normalized_const != 0)
  
  # calculate the approximate (using only the diagonal) of the  inverse of the negative Hessian  
  if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))) {
    temp4  <- .coeff_var(net_stat$node_degree, f, normalized_const,net_stat$m_t,net_stat$offset_tk, net_stat$start_deg + net_stat$g) 
    temp4 <- temp4[non_zero]
  }
  else {
    temp4 <- colSums(net_stat$n_tk[time_non_zero,non_zero_theta, drop = FALSE]^2 * net_stat$m_t[time_non_zero] / normalized_const[time_non_zero]^2)
  }
  if (0 == mode_reg_A) {
    aa <- 1 / (net_stat$sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2 + - temp4 + 
                 lambda * hessian_of_regularization(theta[non_zero_theta]))
    bb <-  1 / (net_stat$sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2 +  
                  lambda * hessian_of_regularization(theta[non_zero_theta])) 
  } else if (1 == mode_reg_A) {
    # mode_reg_A = 1
    
    upper_aa <- net_stat$sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2 + - temp4 
    upper_bb <- net_stat$sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2
   
    if (length(theta) > 0)
      if (length(non_zero_theta) > 0)
        if (length(ok_log10) > 0) {
          upper_aa[ok_log10] <- upper_aa[ok_log10] + lambda * hessian_of_regularization_log10(theta[non_zero_theta[ok_log10]])
          upper_bb[ok_log10] <- upper_bb[ok_log10] + lambda * hessian_of_regularization_log10(theta[non_zero_theta[ok_log10]])
        }
    aa <-  1 / (upper_aa)
    bb <-  1 / (upper_bb)   
  } else if (2 == mode_reg_A) {
    upper_aa <- net_stat$sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2 + - temp4 
    temp_aa <- lambda * hessian_of_regularization_mode_2(theta)
   
    
    if (length(temp_aa) > 0)
      upper_aa <- upper_aa + temp_aa
    
    upper_bb <- net_stat$sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2
    if (length(temp_aa) > 0)
      upper_bb <- upper_bb + lambda * hessian_of_regularization_mode_2(theta)
    aa <-  1 / (upper_aa)
    bb <-  1 / (upper_bb)    
  }
  
  if (length(aa) > 0 && length(bb) > 0)
    cov_bin1[non_zero_theta] <- ifelse(aa > 10^-10, aa, bb) 
 
  cov_bin1[cov_bin1 == Inf] <- 0
  #}  
  
  cov_bin                 <- cov_bin1
  cov_bin[is.na(cov_bin)] <- 0
  cov_bin[cov_bin < 0]    <- 0
  var_log10_bin             <- cov_bin / ifelse(theta != 0 , theta ^ 2 , 1)
  upper_bin               <- exp(log(theta) + 2 * sqrt(var_log10_bin))
  lower_bin               <- exp(log(theta) - 2 * sqrt(var_log10_bin))
  
  
  non_zero_center         <- center_k > 0 & theta > 0 &  var_log10_bin > 0 
  
  
  if (((only_f == FALSE) && (mode_f[1] != "Log_linear")) && (sum(non_zero_center) > 1))  {
    
    
    linear_fit         <- lm(log10(theta[non_zero_center]) ~ log10(center_k[non_zero_center]) , weights = 1 / 
                               var_log10_bin[non_zero_center])
    
    
    
    names(linear_fit$coefficients) <- c("offset","Attachment exponent")
    res                    <- df.residual(linear_fit)
    if (res > 0)
      ci            <- confint(linear_fit,"Attachment exponent")
    else ci         <-"N" 
    
    alpha           <- linear_fit$coefficients[2]
    
  } else if ((only_f == FALSE) && (mode_f[1] != "Log_linear")) {
   
    non_zero_center           <- center_k > 0 & theta > 0 
    linear_fit         <- lm(log10(theta[non_zero_center]) ~ log10(center_k[non_zero_center]))
    
    
    
    names(linear_fit$coefficients) <- c("offset","Attachment exponent")
    res                    <- df.residual(linear_fit)
    if (res > 0)
      ci         <- confint(linear_fit,"Attachment exponent")
    else ci          <-"N" 
    
    alpha           <- linear_fit$coefficients[2]
  } else if (mode_f[1] == "Log_linear") {
    linear_fit <- c(-Inf,-Inf)
    ci          <-"N" 
  } else {
    linear_fit <- c(-Inf,-Inf)
    ci          <-"N" 
    alpha       <- NULL
  } 
  
  
  #non_zero_center <- theta > 0 
 
  #if ((FALSE == only_f) && (mode_f[1] != "log10_linear")) { 
  #     if (sum(non_zero_center) > 0)
  #         alpha_center <- lm(log10(theta[non_zero_center]) ~ log10(center2_k[non_zero_center]))$coefficients[2]
  #     else alpha_center <- NULL   
  # }
  # else 
  #     alpha_center <- NULL   
  ############# Return theta to A #####################################
  A                                  <- rep(0 , net_stat$deg_max)
  cov                                <- rep(0 , net_stat$deg_max)   
  weight_A                           <- rep(1 , net_stat$deg_max)
  A[1 : (net_stat$start_deg + 1)]        <- theta[1:(net_stat$start_deg + 1)] 
  cov[1 : (net_stat$start_deg + 1)]      <- cov_bin[1:(net_stat$start_deg + 1)]
  weight_A[1 : (net_stat$start_deg + 1)] <- 1
  for (i in 1 : net_stat$g) {
    weight_A[(net_stat$begin_deg[i] : net_stat$end_deg[i]) + 1] <- net_stat$interval_length[i]             
    A[(net_stat$begin_deg[i] : net_stat$end_deg[i]) + 1]        <- theta[net_stat$start_deg + i]
    cov[(net_stat$begin_deg[i] : net_stat$end_deg[i])+1]        <- cov_bin[net_stat$start_deg + i] #* weight_A[net_stat$begin_deg[i] + 1]
  }
  interval     <- 0 : (length(A) - 1)
  
  #return(list(k = interval, A = A))
  non_zero     <- which(A > 10^-20)
  non_zero     <- non_zero[non_zero >= net_stat$deg_thresh]
  k_non_zero   <- interval[non_zero]
  
  A            <- A[non_zero] 
  #cc           <- exp(mean(log(k_non_zero[k_non_zero > 0])) - mean(log(A[A > 0])))
  cc           <- 1
  A            <- cc * A
  weight_A     <- weight_A[non_zero]
  cov          <- cc ^ 2 * cov[non_zero]  
  ############### fitting A_k = k^alpha ##################
  var_log10     <- cov / (A ^ 2)
  sd_log10      <- sqrt(var_log10)
  log10_A       <- log(A)
  log10_k       <- log(k_non_zero)
  
  
  ok_var_log10  <- var_log10 > 0
 
  upper_A         <- exp(log(A) + 2 * sd_log10)
  lower_A         <- exp(log(A) - 2 * sd_log10)
  
  
  
  
  ##### Variance of f ####################
  
  #if (rate <= 1 & shape <= 1)
  #   f[-non_zero_f] <- 0
  
  f_new                                        <- rep(offset,net_stat$N)
  names(f_new)                                 <- as.integer(net_stat$node_id)
 
  f_new[as.character(as.integer(net_stat$f_position))]     <- f
  cov_f_new                                    <- rep(0,net_stat$N)
  names(cov_f_new)                             <- as.integer(net_stat$node_id)
  cov_f_new[as.character(as.integer(net_stat$f_position))] <- abs(cov_f)
  non_zero_f                                   <- f_new > 10^-10 & cov_f_new > 10^-10
  
  upper_f                                  <- rep(0,net_stat$N)
  upper_f[non_zero_f]                      <- exp(log(f_new[non_zero_f]) + 2 * sqrt(cov_f_new[non_zero_f] / f_new[non_zero_f] ^ 2))
  
  lower_f                                  <- rep(0,net_stat$N)
  lower_f[non_zero_f]                      <- exp(log(f_new[non_zero_f]) - 2 * sqrt(cov_f_new[non_zero_f] / f_new[non_zero_f] ^ 2))
  
  if (FALSE == only_f) {
    if (mode_f[1] != "Log_linear")
      alpha = linear_fit$coefficients[2]
    else names(alpha) <- "Estimated attachment exponent"
  }
  
  
  
  
  
  result <- list(# estimated PA function and its variances, confidence interval 
    k       = k_non_zero ,  A             = A          , var_A       = cov       , var_logA = var_log10 ,
    upper_A = upper_A    ,  lower_A       = lower_A    , center_k = center_k,  
    theta   = theta      ,  upper_bin     = upper_bin  , lower_bin   = lower_bin , var_bin  = cov_bin ,var_logbin = var_log10_bin,
    
    # estimated attachment exponent alpha, and the log-linear fit 
    alpha   = alpha      ,  loglinear_fit = linear_fit , ci          = ci        ,
    
    
    alpha_series = ifelse(rep(mode_f[1] == "Log_linear", length(alpha_series)),alpha_series,-1),
    
    # estimated node fitnesses and their variances, confidence intervals
    f        = f_new     ,  var_f         = cov_f_new  , 
    upper_f  = upper_f   ,  lower_f       = lower_f    , # confidence intervals
    
    # values of the objective function over iterations
    objective_value = log10_likelihood,
    
    # other parameters specified
    mode_f = mode_f[1] , true_A = true_A , true_f = true_f,
    only_PA = only_PA  , only_f = only_f , lambda = lambda, 
    shape = shape, rate = rate, 
    deg_threshold = net_stat$deg_thresh, stop_cond = stop_cond, auto_lambda = auto_lambda, ratio = ratio, 
    g = net_stat$g, diverge_zero = diverge_zero
  )
  class(result) <- "PAFit_result"
  return(result)
}
