# function to estimate jointly the attachment function and node fitness  
PAFit <- function(net_stat, 
                  only_PA        = FALSE       , only_f         = FALSE        , 
                  mode_f         = "Linear_PA" ,
                  true_A         = NULL        , true_f         = NULL         , 
                  
                  mode_reg_A     = 0           , weight_PA_mode = 1            ,
                  s              = 10          , lambda         = 1            , 
                  auto_lambda    = TRUE        , r              = 0.01         ,

                  alpha_start    = 1           , start_mode_A   = "Log_linear" , 
                  start_mode_f   = "Constant"  ,
                  
                  
                  auto_stop      = TRUE        , stop_cond      = 10^-7        , 
                  iteration      = 200         , max_iter       = 200000       , 
                  debug          = FALSE       , q              = 1            , 
                  step_size      = 0.5         ,
                  
                  normalized_f   = FALSE       , interpolate = TRUE) {
    if ((net_stat$only_PA == TRUE) & (only_PA == FALSE)) {
        warning("The net_stat does not support estimation of node fitness. It will be run with option 'only_PA = TRUE'.
                 Please re-run GetStatistics again with the option 'only_PA = FALSE' if you also want to estimate fitnesses.")
        only_PA <- TRUE  
    }

    shape          <- s
    rate           <- s
    ratio          <- r
    if (s <= 0)
        stop("Error: shape and rate should be positive number")
  
    non_zero_theta     <- which(net_stat$Sum_m_k > 0)
    num_nonzero        <- length(non_zero_theta)
    theta              <- rep(1,length(net_stat$Sum_m_k))
    
    
    
     for (ii in 1:length(theta))
        if ((mode_f[1] == "Constant_PA") & (only_f == TRUE)) 
        {
            theta[ii] <- 1;  
        } else if ((only_f == TRUE) & (!is.null(true_A))) {
            theta[ii] <- true_A[ii];    
        }
        else if ("Log_linear" == start_mode_A[1]){    
            if (ii <= net_stat$start_deg & ii %in% non_zero_theta)
               theta[ii] <- ii^alpha_start
            else 
            if (ii > net_stat$start_deg & ii %in% non_zero_theta)
               theta[ii] <- net_stat$begin_deg[ii - net_stat$start_deg]^alpha_start
        } else if ("Random" == start_mode_A[1]) {
            theta <- runif(length(theta),0,1);  
        }
    
    alpha <- 1.0
    if (theta[length(theta)] == 1 || theta[length(theta)] == 0)
        theta[length(theta)] <- theta[length(theta) - 1] + 1
    #print(theta)
    theta[which(theta <= 0)] <- 1
    if (mode_f[1] != "Log_linear")
        theta[-non_zero_theta] <- 0
    if (mode_f[1] == "Log_linear") {
      alpha <- 1.0 * alpha_start; #starting value of alpha
    }
    else theta         <- theta/sum(theta)
    #print(theta)
    #if (include_zero == 0)
    non_zero_f    <- which(net_stat$z_j >= 0)

    #if (shape <= 1 & rate <= 1)
    #    non_zero_f <- which(net_stat$z_j > 0)
    if ("Constant" == start_mode_f)
        f             <- rep(2,length(net_stat$f_position))
    else if ("Random" == start_mode_f)
        f             <- rgamma(n = length(net_stat$f_position), shape = s, rate = s);
    f                 <- length(f) * f / sum(f)
    #if (shape < 1 & rate < 1)
    #     f[net_stat$z_j == 0] <- 0
    
    if (TRUE == only_PA) {
        if (is.null(true_f))
            f[] <- 1
        else
            f[] <- true_f  
    }
    log_likelihood   <- vector()
    # get the center of each bin
    center_k    <- rep(0, length(theta)) 
    center2_k   <- center_k
    if (net_stat$start_deg > 0) {
      center_k[1:net_stat$start_deg]  <- 0:(net_stat$start_deg - 1)
      center2_k[1:net_stat$start_deg] <- 0:(net_stat$start_deg - 1)
    }
    for (i in 1:net_stat$G) {
      if (net_stat$begin_deg[i] != 0) {
        # center_k[i]  <- round((net_stat$begin_deg[i] + net_stat$end_deg[i])/2)  
        #center_k[net_stat$start_deg + i] <- round(net_stat$begin_deg[i]*sqrt((net_stat$begin_deg[i] + net_stat$interval_length[i] - 1)/ net_stat$begin_deg[i]))
        #center_k[net_stat$start_deg + i] <- net_stat$begin_deg[i]
        center_k[net_stat$start_deg + i] <- net_stat$end_deg[i]
        #center2_k[net_stat$start_deg + i] <- round((net_stat$begin_deg[i] + net_stat$end_deg[i])/2)  
        center2_k[net_stat$start_deg + i] <- net_stat$end_deg[i]
      } else {
        center_k[net_stat$start_deg + i] <- net_stat$end_deg[i]  
        #center2_k[i]  <- round((net_stat$begin_deg[i] + net_stat$end_deg[i])/2)   
        center2_k[net_stat$start_deg + i] <- net_stat$end_deg[i]
      }
      #        
    }
    #print(center_k)
      

    if (0 == mode_reg_A) {
        update_theta   <- non_zero_theta[-c(1 , 2 , num_nonzero , num_nonzero - 1)]
        noupdate_theta <- non_zero_theta[c(1 , 2 , num_nonzero , num_nonzero - 1)]
        minus_1        <- non_zero_theta[-c(1 , num_nonzero , num_nonzero - 1 , num_nonzero - 2)]
        minus_2        <- non_zero_theta[-c(num_nonzero , num_nonzero - 1 , num_nonzero - 2 , num_nonzero - 3)] 
        plus_1         <- non_zero_theta[-c(1 , 2 , 3 , num_nonzero)]
        plus_2         <- non_zero_theta[-c(1 , 2 , 3 , 4)] 
        #num_ok_log     <- 1
    } else if (1 == mode_reg_A) {
          ok_log     <- which(center_k[non_zero_theta] > 1)
          not_ok_log <- which(center_k[non_zero_theta] <= 1)
          num_ok_log <- length(ok_log) 
          #print(not_ok_log)
          update_theta   <- non_zero_theta[ok_log[-c(1, 2, num_ok_log, num_ok_log - 1)]]
          noupdate_theta <- non_zero_theta[ok_log[c(1, 2, num_ok_log, num_ok_log - 1)]]
          minus_1 <- non_zero_theta[ok_log[-c(1, num_ok_log, num_ok_log - 1, num_ok_log - 2)]]
          minus_2 <- non_zero_theta[ok_log[-c(num_ok_log, num_ok_log - 1, num_ok_log - 2, num_ok_log - 3)]] 
          plus_1  <- non_zero_theta[ok_log[-c(1, 2, 3, num_ok_log)]]
          plus_2  <- non_zero_theta[ok_log[-c(1, 2, 3, 4)]] 
    } else if (2 == mode_reg_A) {
          ok_log     <- which(center_k[non_zero_theta] > 0)
          not_ok_log <- which(center_k[non_zero_theta] == 0)
          num_ok_log <- length(ok_log) 
          #print(ok_log)
          #print(not_ok_log)
          #print(non_zero_theta)
          #print(center_k)
          if (num_ok_log >= 3) { 
              update_theta   <- non_zero_theta[ok_log[-c(1, 2, num_ok_log, num_ok_log - 1)]]
              noupdate_theta <- non_zero_theta[ok_log[c(1, 2, num_ok_log, num_ok_log - 1)]]
              minus_1 <- non_zero_theta[ok_log[-c(1, num_ok_log, num_ok_log - 1, num_ok_log - 2)]]
              minus_2 <- non_zero_theta[ok_log[-c(num_ok_log, num_ok_log - 1, num_ok_log - 2, num_ok_log - 3)]] 
              plus_1  <- non_zero_theta[ok_log[-c(1, 2, 3, num_ok_log)]]
              plus_2  <- non_zero_theta[ok_log[-c(1, 2, 3, 4)]] 
              u       <- matrix(0,nrow = num_ok_log, ncol = num_ok_log)
          #print(update_theta)
          #print(minus_1)
          
              for (k1 in 2:(num_ok_log - 1)) {
                  u[k1,k1]       <- 1/(log(center_k[non_zero_theta[ok_log[k1 + 1 ]]]) - log(center_k[non_zero_theta[ok_log[k1]]])) +
                                            1/(log(center_k[non_zero_theta[ok_log[k1]]]) - log(center_k[non_zero_theta[ok_log[k1 - 1]]])) 
                  u[k1,k1 - 1]   <- 1/(log(center_k[non_zero_theta[ok_log[k1]]]) - log(center_k[non_zero_theta[ok_log[k1 - 1]]]))
                  u[k1,k1 + 1]   <- 1/(log(center_k[non_zero_theta[ok_log[k1 + 1 ]]]) - log(center_k[non_zero_theta[ok_log[k1]]]))
              }
          #print(u)
          extract_u <- function(i,j) {
            result <- NULL
            for (index in 1:length(i))
              result <- c(result,u[i[index],j[index]])  
            return(result)
          }
          update_k <- 3:(num_ok_log-2)
         }
          #print(center_k[non_zero_theta[ok_log]])
          #print(u)     
    }
    
    PA_offset <- 1
    if(theta[1] != 0)
        theta[1]  <- PA_offset
    
    #starting value for the offset
    offset <- 1
    #weights of the regularization term
    if (weight_PA_mode[1] == 0)
        w_k <- net_stat$Sum_m_k/sum(net_stat$Sum_m_k)
    else {
        if (weight_PA_mode[1] == 1)
            w_k <- rep(1/length(theta),length(theta))
        else if (weight_PA_mode[1] == 2)
            w_k <- 1/net_stat$Sum_m_k/sum(1/net_stat$Sum_m_k)
    }
    
    if (TRUE == auto_lambda) {
       #lambda <- ratio   
        lambda <- ratio * sum(net_stat$Sum_m_k)
    }
    if (TRUE == auto_stop)
        iteration <- max_iter

    normalized_const <- rep(0, dim(net_stat$n_tk)[1])
    alpha_series <- c()
    
    break_flag  <- FALSE
    count_break <- 0
    
    cal_reg_A <- function() {
    #print(theta)  
      if (lambda > 0) {  
        if (0 == mode_reg_A) {
            return(sum(lambda * w_k[non_zero_theta[-c(1,num_nonzero)]] * (log(theta[non_zero_theta[-c(1 , 2)]]) + 
                                                                          log(theta[non_zero_theta[-c(num_nonzero , num_nonzero - 1)]]) - 
                                                                          2 * log(theta[non_zero_theta[-c(1 , num_nonzero)]]))^2))
        } else if (1 == mode_reg_A) {
              return(sum(lambda * w_k[non_zero_theta[ok_log[-c(1 , num_ok_log)]]] *
                       (log(theta[non_zero_theta[ok_log[-c(1 , 2)]]]) / log(center_k[non_zero_theta[ok_log[-c(1 , 2)]]]) + 
                          log(theta[non_zero_theta[ok_log[-c(num_ok_log , num_ok_log - 1)]]]) / 
                          log(center_k[non_zero_theta[ok_log[-c(num_ok_log , num_ok_log - 1)]]])  - 
                          2 * log(theta[non_zero_theta[ok_log[-c(1 , num_ok_log)]]]) / 
                          log(center_k[non_zero_theta[ok_log[-c(1 , num_ok_log)]]]) ) ^ 2))
        } else if (2 == mode_reg_A) {
            #print(w_k[non_zero_theta[ok_log[-c(1,num_ok_log)]]])  
            #print(log(theta[non_zero_theta[ok_log[-c(1,2)]]]))
            #print((
            #  (log(theta[non_zero_theta[ok_log[-c(1,2)]]]) -  log(theta[non_zero_theta[ok_log[-c(1,num_ok_log)]]])) / 
            #    (log(center_k[non_zero_theta[ok_log[-c(1,2)]]]) -  log(center_k[non_zero_theta[ok_log[-c(1,num_ok_log)]]]))
            #))
            #print((
            #  (log(theta[non_zero_theta[ok_log[-c(1,num_ok_log)]]]) -  log(theta[non_zero_theta[ok_log[-c(num_ok_log,num_ok_log - 1)]]])) / 
            #    (log(center_k[non_zero_theta[ok_log[-c(1,num_ok_log)]]]) -  log(center_k[non_zero_theta[ok_log[-c(num_ok_log, num_ok_log - 1)]]])) 
            #) )
            #print(   w_k[non_zero_theta[not_ok_log]] * 
            #        (log(theta[non_zero_theta[not_ok_log]]) - log(theta[non_zero_theta[ok_log[1]]])) )
          if (num_ok_log >= 3) {
              return_temp <- sum(lambda * w_k[non_zero_theta[ok_log[-c(1,num_ok_log)]]] *
                               ( (
                                 (log(theta[non_zero_theta[ok_log[-c(1,2)]]]) -  log(theta[non_zero_theta[ok_log[-c(1,num_ok_log)]]])) / 
                                   (log(center_k[non_zero_theta[ok_log[-c(1,2)]]]) -  log(center_k[non_zero_theta[ok_log[-c(1,num_ok_log)]]]))
                               ) -  
                                 (
                                   (log(theta[non_zero_theta[ok_log[-c(1,num_ok_log)]]]) -  log(theta[non_zero_theta[ok_log[-c(num_ok_log,num_ok_log - 1)]]])) / 
                                     (log(center_k[non_zero_theta[ok_log[-c(1,num_ok_log)]]]) -  log(center_k[non_zero_theta[ok_log[-c(num_ok_log, num_ok_log - 1)]]])) 
                                 ) 
                               )^2)
              if (length(not_ok_log) > 0)
                  return_temp <- return_temp +  lambda * w_k[non_zero_theta[not_ok_log]] * 
                                 (log(theta[non_zero_theta[not_ok_log]]) - log(theta[non_zero_theta[ok_log[1]]])) ^2    
              return(return_temp);  
          } else return(0);
          
        }
      } else return(0);
    }
    
    objective_function_value <- function(theta,f,offset,normalized_const) {
        time_non_zero       <- which(normalized_const > 0)
        non_zero_f_temp     <- which(f > 0)
        non_zero_theta_temp <- which(theta > 0) 
        #print(time_non_zero)
        #print(non_zero_theta_temp)
        if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))) {
          value <- sum(net_stat$z_j[non_zero_f_temp] * log(f[non_zero_f_temp])) + 
                                 sum(net_stat$offset_m_tk)*log(offset) +  
                                 (shape - 1) * log(offset) - rate * offset +    
                                 sum(net_stat$Sum_m_k[non_zero_theta_temp] * log(theta[non_zero_theta_temp])) - 
                                 sum(net_stat$m_t[time_non_zero] * log(normalized_const[time_non_zero])) + 
                                 (shape - 1) * (sum(log(f[non_zero_f_temp]))) - rate * sum(f[non_zero_f_temp]) - cal_reg_A()
        #print(log_likelihood)
      }
      else if (is.null(true_f) && (TRUE == only_PA))  { 
          #print("Go inside here!")  
          value <- sum(net_stat$Sum_m_k[non_zero_theta_temp] * log(theta[non_zero_theta_temp])) - 
                               sum(net_stat$m_t[time_non_zero] * log(normalized_const[time_non_zero]))
          
          #print(paste0("value: ",value));  
          value <- value - cal_reg_A()
      }
      else  if ((is.null(true_A) && (TRUE == only_f)))
          value <- sum(net_stat$z_j[non_zero_f_temp] * log(f[non_zero_f_temp])) - 
                                 sum(net_stat$m_t[time_non_zero] * log(normalized_const[time_non_zero])) + 
                                 sum(net_stat$offset_m_tk)*log(offset) + 
                                 (shape - 1) * log(offset) - rate * offset +    
                                 (shape - 1) * (sum(log(f[non_zero_f_temp]))) - rate * sum(f[non_zero_f_temp])  
       
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
        #print(normalized_const)  
        #print(cal_reg_A());  
        #print(normalized_const);
        #print(theta);
        #print(log(theta[non_zero_theta]));
        #print(offset)
        
        log_likelihood    <- c(log_likelihood,objective_function_value(theta,f,offset,normalized_const))  ;
        
        time_non_zero     <- which(normalized_const != 0)
        non_zero_f_temp   <- which(f > 0)
        

        if (TRUE == debug) {
            print(log_likelihood[length(log_likelihood)])
        }
        if (length(log_likelihood) > 1)
            if (log_likelihood[length(log_likelihood)] < log_likelihood[length(log_likelihood) - 1])
                stop("Warning: Log likelihood decreased.")  
        
        
        if (TRUE == auto_stop)
            if (length(log_likelihood) > 1)
                tryCatch({if (abs(log_likelihood[length(log_likelihood)] - log_likelihood[length(log_likelihood) - 1]) / 
                          (abs(log_likelihood[length(log_likelihood) - 1]) + 1) < stop_cond)
                           break_flag <- TRUE;},error = function(e) { #print(as.vector(normalized_const));print(f[non_zero_f]);
                                                                      #print(non_zero_f);
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
                                          error = function(e) {#print("Problem with inversion"); 
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
                log_candidate <- objective_function_value(theta_temp,f_temp,offset_temp,normalized_const_temp)
                if (log_candidate > log_likelihood[length(log_likelihood)]) {
                    theta  <- candidate[1 : length(theta)]
                    f      <- candidate[(length(theta) + 1) : (length(theta) + length(f))]
                    offset <- candidate[length(candidate)]  
                    log_likelihood[length(log_likelihood)]  <- log_candidate
                    normalized_const                        <- normalized_const_temp
                }
                else {
                #print(paste("ll of rejected candidate:",log_candidate))
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
        }
        }
        ############### End of quasi-Newton acceleration ################
        
        
        if (break_flag) {
            break;   
        } else { 
          break_flag  <- FALSE
        }

        #print(offset)
        if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))){
          .normalized_constant(normalized_const,net_stat$node_degree,theta,f,net_stat$offset_tk,offset) 
        }
        else normalized_const <- as.vector(net_stat$n_tk[,non_zero_theta]%*%theta[non_zero_theta])
        #.normalized_constant(normalized_const, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
        
        #print(normalized_const)
        
       
       
        ##################### Update f ######################
        if (FALSE == only_PA){
            
            .update_f(f,non_zero_f,net_stat$node_degree,theta,net_stat$z_j,normalized_const,net_stat$m_t,shape,rate,offset)
            
          #print(f)
          
          ### NORMALIZED f
          #normalize_f <- mean(c(f,offset))
          #f           <- f/normalize_f
          #offset      <- offset/normalize_f
          
          
          
            #update offset
            #if (net_stat$deg_thresh > 0)
            .normalized_constant(normalized_const, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
            offset <- .update_offset(net_stat$offset_tk, net_stat$offset_m_tk, theta, normalized_const,net_stat$m_t, shape,rate)
            #if (debug_PA == TRUE)
            #    print(paste0("Offset: ",offset))
            .normalized_constant(normalized_const, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
            time_non_zero     <- which(normalized_const != 0)
            #print(f)
        }

        #### Update PA_offset if mode_f = "Linear_PA" #####
        #if ((TRUE == only_f) && (mode_f[1] == "Linear_PA")) {
        #    PA_offset <- .update_PA_offset(normalized_const,f,net_stat$node_degree,net_stat$m_t,net_stat$Sum_m_k,
        #                                   net_stat$offset_tk);
        #    theta[1]  <- PA_offset
        #}
        #####################  Update A #######################
       
        if (FALSE == only_f) {
            if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))) {
                #print("problem in coeff theta")
                temp5  <- .coeff_theta(net_stat$node_degree, f, normalized_const,net_stat$m_t,net_stat$start_deg + net_stat$G)   
                if (length(time_non_zero) > 1) {
                    temp4  <- temp5[non_zero_theta] + colSums(net_stat$m_t[time_non_zero] / normalized_const[time_non_zero] * offset *
                                                          net_stat$offset_tk[time_non_zero,non_zero_theta])
                } else {
                    temp4  <- temp5[non_zero_theta] + net_stat$m_t[time_non_zero] / normalized_const[time_non_zero] * offset *
                                                              net_stat$offset_tk[time_non_zero,non_zero_theta]  
                }
            }
            else {
                if (length(time_non_zero) > 1) {  
                    temp4 <- colSums(net_stat$m_t[time_non_zero]/normalized_const[time_non_zero] * 
                                     net_stat$n_tk[time_non_zero,non_zero_theta])
                } else {
                    temp4 <- net_stat$m_t[time_non_zero]/normalized_const[time_non_zero] * 
                                     net_stat$n_tk[time_non_zero,non_zero_theta]  
                }
            }
           #print(temp4)
           #print(ok_log)
           #print("reach here")
            if ((lambda <= 0) || ((0 != mode_reg_A) && (num_ok_log < 3))) {
              
                    theta[non_zero_theta] <- net_stat$Sum_m_k[non_zero_theta] / temp4
            }
            else if (0 == mode_reg_A) {
                g_1  <- function(x) {
                        (net_stat$Sum_m_k[non_zero_theta[1]] - 2*w_k[non_zero_theta[2]] * lambda * (-2 * log(theta[non_zero_theta[2]]) + 
                        log(theta[non_zero_theta[3]]) - 3 * log(theta[non_zero_theta[1]])))/x -   
                         temp4[1] - 8 * lambda * w_k[non_zero_theta[2]] * log(x) / x}
                
                g_2 <- function(x) {
                       (net_stat$Sum_m_k[non_zero_theta[2]] - 2*lambda * 
                          (2*w_k[non_zero_theta[2]]*(-2*log(theta[non_zero_theta[2]]) -log(theta[non_zero_theta[3]]) - log(theta[non_zero_theta[1]])) + 
                             w_k[non_zero_theta[3]]*(-2*log(theta[non_zero_theta[3]]) +log(theta[non_zero_theta[4]]) -3*log(theta[non_zero_theta[2]]) ))) / x - 
                         temp4[2] - (16*w_k[non_zero_theta[2]] + 8*w_k[non_zero_theta[3]]) * lambda * log(x) / x}
               g_semiend <- function(x) {
                            (net_stat$Sum_m_k[non_zero_theta[num_nonzero - 1]] - 2*lambda * 
                            (2*w_k[non_zero_theta[num_nonzero - 1]]*(-2*log(theta[non_zero_theta[num_nonzero - 1]]) - log(theta[non_zero_theta[num_nonzero]]) - log(theta[non_zero_theta[num_nonzero - 2]])) +
                             w_k[non_zero_theta[num_nonzero -2]] *(-3*log(theta[non_zero_theta[num_nonzero - 1]]) + log(theta[non_zero_theta[num_nonzero - 3]]) - 2*log(theta[non_zero_theta[num_nonzero - 2]])))) / x - 
                              temp4[num_nonzero - 1] - 
                              (16*w_k[non_zero_theta[num_nonzero - 1]] + 8*w_k[non_zero_theta[num_nonzero -2]]) * lambda * log(x) / x}
               g_end     <-  function(x) {
                             (net_stat$Sum_m_k[non_zero_theta[num_nonzero]] - 2*lambda * w_k[non_zero_theta[num_nonzero -1]] *(-3*log(theta[non_zero_theta[num_nonzero]]) + log(theta[non_zero_theta[num_nonzero - 2]]) - 2*log(theta[non_zero_theta[num_nonzero - 1]]))) / x  -   
                             temp4[num_nonzero] - 8*w_k[non_zero_theta[num_nonzero -1]]*lambda*log(x)/x}  

               U_1 <- net_stat$Sum_m_k[update_theta] - 2 * lambda * (2*w_k[update_theta]*(-2*log(theta[update_theta]) - log(theta[plus_1]) - log(theta[minus_1])) + 
                                                                   w_k[minus_1]*(-3*log(theta[update_theta]) + log(theta[minus_2]) - 2*log(theta[minus_1])) + 
                                                                   w_k[plus_1]*(-2*log(theta[plus_1]) + log(theta[plus_2]) - 3*log(theta[update_theta])))
               U_2 <- (16*w_k[update_theta] + 8*w_k[minus_1] + 8*w_k[plus_1]) * lambda
               U_3 <- temp4[-c(1,2,num_nonzero-1,num_nonzero)]
               theta[non_zero_theta[1]] <- tryCatch(uniroot(g_1,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                    error = function(e) {return(theta[non_zero_theta[1]])})
        
               theta[non_zero_theta[2]] <- tryCatch(uniroot(g_2,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                   error = function(e) {return(theta[non_zero_theta[2]])})
               #parallelization here
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
                (net_stat$Sum_m_k[non_zero_theta[ok_log[1]]] - 2 * w_k[non_zero_theta[ok_log[2]]] * lambda * (-2 * log(theta[non_zero_theta[ok_log[2]]]) / 
                                                                                                          log(center_k[non_zero_theta[ok_log[2]]])  + 
                                                                                          log(theta[non_zero_theta[ok_log[3]]]) / log(center_k[non_zero_theta[ok_log[3]]]) - 
                                                                                          3 * log(theta[non_zero_theta[ok_log[1]]]) / log(center_k[non_zero_theta[ok_log[1]]])  ))/x -   
                  temp4[ok_log[1]] - 8 * lambda * w_k[non_zero_theta[ok_log[2]]] * log(x) / x / log(center_k[non_zero_theta[ok_log[1]]])}
              
              g_2 <- function(x) {
                (net_stat$Sum_m_k[non_zero_theta[ok_log[2]]] - 2 * lambda * 
                   (2 * w_k[non_zero_theta[ok_log[2]]] * (- 2 * log(theta[non_zero_theta[ok_log[2]]]) / log(center_k[non_zero_theta[ok_log[2]]]) - 
                                                            log(theta[non_zero_theta[ok_log[3]]]) / log(center_k[non_zero_theta[ok_log[3]]]) - 
                                                            log(theta[non_zero_theta[ok_log[1]]]) / log(center_k[non_zero_theta[ok_log[1]]])   ) + 
                      w_k[non_zero_theta[ok_log[3]]]*(- 2 * log(theta[non_zero_theta[ok_log[3]]]) / log(center_k[non_zero_theta[ok_log[3]]])  +
                                                        log(theta[non_zero_theta[ok_log[4]]]) / log(center_k[non_zero_theta[ok_log[4]]]) -
                                                        3 * log(theta[non_zero_theta[ok_log[2]]]) / log(center_k[non_zero_theta[ok_log[2]]])  ))) / x - 
                  temp4[ok_log[2]] - (16*w_k[non_zero_theta[ok_log[2]]] + 8*w_k[non_zero_theta[ok_log[3]]]) * lambda * log(x) / x / log(center_k[non_zero_theta[ok_log[2]]])}
          
                  g_semiend <- function(x) {
                (net_stat$Sum_m_k[non_zero_theta[ok_log[num_ok_log - 1]]] - 2*lambda * 
                   (2*w_k[non_zero_theta[ok_log[num_ok_log - 1]]] * ( -2 * log(theta[non_zero_theta[ok_log[num_ok_log - 1]]]) / log(center_k[non_zero_theta[ok_log[num_ok_log - 1]]]) - 
                                                                        log(theta[non_zero_theta[ok_log[num_ok_log]]]) / log(center_k[non_zero_theta[ok_log[num_ok_log]]]) - 
                                                                  log(theta[non_zero_theta[ok_log[num_ok_log - 2]]]) / log(center_k[non_zero_theta[ok_log[num_ok_log - 2]]]) ) +
                      w_k[non_zero_theta[ok_log[num_ok_log -2]]] *( - 3 * log(theta[non_zero_theta[ok_log[num_ok_log - 1]]]) / log(center_k[non_zero_theta[ok_log[num_ok_log - 1]]]) + 
                                                                      log(theta[non_zero_theta[ok_log[num_ok_log - 3]]]) / log(center_k[non_zero_theta[ok_log[num_ok_log - 3]]]) - 
                                                                      2 * log(theta[non_zero_theta[ok_log[num_ok_log - 2]]]) / log(center_k[non_zero_theta[ok_log[num_ok_log - 2]]]) ) ) ) / x - 
                  temp4[ok_log[num_ok_log - 1]] - 
                  (16*w_k[non_zero_theta[ok_log[num_ok_log - 1]]] + 8 * w_k[non_zero_theta[ok_log[num_ok_log -2]]]) * lambda * log(x) / x / log(center_k[non_zero_theta[ok_log[num_ok_log - 1]]])}
              
                  g_end     <-  function(x) {
                (net_stat$Sum_m_k[non_zero_theta[ok_log[num_ok_log]]] - 2 * lambda * w_k[non_zero_theta[ok_log[num_ok_log -1]]] * (-3 * 
                                               log(theta[non_zero_theta[ok_log[num_ok_log]]]) / log(center_k[non_zero_theta[ok_log[num_ok_log]]]) + 
                                               log(theta[non_zero_theta[ok_log[num_ok_log - 2]]]) / log(center_k[non_zero_theta[ok_log[num_ok_log - 2]]]) - 
                                               2 * log(theta[non_zero_theta[ok_log[num_ok_log - 1]]]) / log(center_k[non_zero_theta[ok_log[num_ok_log - 1]]])   )) / x  -   
                  temp4[ok_log[num_ok_log]] - 8 * w_k[non_zero_theta[ok_log[num_ok_log - 1]]] * lambda * log(x)/x / log(center_k[non_zero_theta[ok_log[num_ok_log]]]) }  
              
            
                  
                  U_1 <- net_stat$Sum_m_k[update_theta] - 2 * lambda * (2 * w_k[update_theta] * (-2 * log(theta[update_theta]) / log(center_k[update_theta]) - 
                                                                                               log(theta[plus_1]) /  log(center_k[plus_1]) - log(theta[minus_1]) / log(center_k[minus_1]) ) + 
                                                                      w_k[minus_1]*(-3 * log(theta[update_theta]) / log(center_k[update_theta]) + 
                                                                                      log(theta[minus_2]) / log(center_k[minus_2]) - 2 * log(theta[minus_1]) / log(center_k[minus_1])  ) + 
                                                                      w_k[plus_1]*(-2 * log(theta[plus_1]) / log(center_k[plus_1]) + log(theta[plus_2]) / log(center_k[plus_2]) - 
                                                                                     3 * log(theta[update_theta]) / log(center_k[update_theta])))
                  
                  U_2 <- (16*w_k[update_theta] + 8*w_k[minus_1] + 8*w_k[plus_1]) * lambda
                  
                  U_3 <- temp4[ok_log[-c(1,2,num_ok_log - 1,num_ok_log)]]
                  
                  
              theta[non_zero_theta[ok_log[1]]] <- tryCatch(uniroot(g_1,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                   error = function(e) {return(theta[non_zero_theta[ok_log[1]]])})
              
              theta[non_zero_theta[ok_log[2]]] <- tryCatch(uniroot(g_2,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                   error = function(e) {return(theta[non_zero_theta[ok_log[2]]])})
              #parallelization here
              for (jj in 1:length(update_theta)) {
                g <- function(x){U_1[jj]/x - U_2[jj]*log(x)/x / log(center_k[update_theta[jj]]) - U_3[jj]}
                theta[update_theta[jj]] <- tryCatch(uniroot(g,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                    error = function(e) theta[update_theta[jj]])
              }
              theta[non_zero_theta[ok_log[num_ok_log - 1]]] <- tryCatch(uniroot(g_semiend,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                                 error = function(e) 
                                                                 {return(theta[non_zero_theta[ok_log[num_ok_log - 1]]])})
              theta[non_zero_theta[ok_log[num_ok_log]]] <- tryCatch(uniroot(g_end,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                             error = function(e) {return(theta[non_zero_theta[ok_log[num_ok_log]]])})
              
              if (length(not_ok_log) > 0) 
                theta[non_zero_theta[not_ok_log]] <- net_stat$Sum_m_k[non_zero_theta[not_ok_log]]/temp4[not_ok_log]  
            } else if (2 == mode_reg_A) {
              #mode_reg_A == 2  

              g_1_new <- function(x) {
                (net_stat$Sum_m_k[non_zero_theta[ok_log[1]]] - 2 * u[2,1] * w_k[non_zero_theta[ok_log[2]]] * lambda * 
                   (- u[2,2] * log(theta[non_zero_theta[ok_log[2]]])   + 
                      u[2,3] * log(theta[non_zero_theta[ok_log[3]]])  - 
                      3 * u[2,1] * log(theta[non_zero_theta[ok_log[1]]])  ) + 
                      2 * lambda * w_k[non_zero_theta[not_ok_log]] * 
                      (log(theta[non_zero_theta[ok_log[1]]])  + log(theta[non_zero_theta[not_ok_log]]) )) / x -   
                  temp4[ok_log[1]] - (8 * lambda * w_k[non_zero_theta[ok_log[2]]] * u[2,1]^2 + 
                                      4 * lambda * w_k[non_zero_theta[not_ok_log]])* log(x) / x }
              
              g_2 <- function(x) {
                (net_stat$Sum_m_k[non_zero_theta[ok_log[2]]] - lambda * 
                   (-2 * u[2 , 2] * w_k[non_zero_theta[ok_log[2]]] * ( u[2,2] * log(theta[non_zero_theta[ok_log[2]]]) + 
                                                                       u[2,3] * log(theta[non_zero_theta[ok_log[3]]]) +  
                                                                       u[2,1] * log(theta[non_zero_theta[ok_log[1]]]) ) +
                      
                      2 * u[3,2] * w_k[non_zero_theta[ok_log[3]]]*(- u[3,3] * log(theta[non_zero_theta[ok_log[3]]])  +
                                                                   u[3, 4] * log(theta[non_zero_theta[ok_log[4]]]) - 
                                                                  3 * u[3,2] * log(theta[non_zero_theta[ok_log[2]]])  ) )) / x - 
                  
                  temp4[ok_log[2]] + (-4 * u[2,2]^2 * w_k[non_zero_theta[ok_log[2]]] - 
                                      8 * u[3,2]^2 * w_k[non_zero_theta[ok_log[3]]]) * lambda * log(x) / x     }
              
              g_semiend <- function(x) {
                (net_stat$Sum_m_k[non_zero_theta[ok_log[num_ok_log - 1]]] - lambda * 
                   (-2 * u[num_ok_log - 1, num_ok_log - 1] * w_k[non_zero_theta[ok_log[num_ok_log - 1]]] * 
                      ( u[num_ok_log - 1, num_ok_log - 1] * log(theta[non_zero_theta[ok_log[num_ok_log - 1]]]) + 
                          u[num_ok_log - 1, num_ok_log] * log(theta[non_zero_theta[ok_log[num_ok_log]]])  + 
                          u[num_ok_log - 1, num_ok_log -  2] * log(theta[non_zero_theta[ok_log[num_ok_log - 2]]]) )  +
                     2 * u[num_ok_log - 2, num_ok_log - 1] * w_k[non_zero_theta[ok_log[num_ok_log -2]]] * 
                      ( - 3 * u[num_ok_log - 2, num_ok_log - 1]* log(theta[non_zero_theta[ok_log[num_ok_log - 1]]]) + 
                          u[num_ok_log - 2, num_ok_log - 3] * log(theta[non_zero_theta[ok_log[num_ok_log - 3]]])  - 
                          u[num_ok_log - 2, num_ok_log - 2] * log(theta[non_zero_theta[ok_log[num_ok_log - 2]]])  ) ) ) / x - 
                  temp4[ok_log[num_ok_log - 1]] - 
                  (4 * u[num_ok_log - 1, num_ok_log - 1]^2 * w_k[non_zero_theta[ok_log[num_ok_log - 1]]] + 8 * u[num_ok_log - 2, num_ok_log - 1] ^ 2 * 
                     w_k[non_zero_theta[ok_log[num_ok_log -2]]]) * lambda * log(x) / x }
              
              g_end     <-  function(x) {
                (net_stat$Sum_m_k[non_zero_theta[ok_log[num_ok_log]]] - 
                   2 * u [num_ok_log - 1, num_ok_log] * lambda * w_k[non_zero_theta[ok_log[num_ok_log -1]]] * 
                   (-3 * u [num_ok_log - 1, num_ok_log] * log(theta[non_zero_theta[ok_log[num_ok_log]]])  + 
                      u [num_ok_log - 1, num_ok_log - 2 ] * log(theta[non_zero_theta[ok_log[num_ok_log - 2]]])  - 
                      u [num_ok_log - 1, num_ok_log - 1] * log(theta[non_zero_theta[ok_log[num_ok_log - 1]]])   )) / x  -   
                  temp4[ok_log[num_ok_log]] - 
                  8 * u [num_ok_log - 1, num_ok_log]^2 * w_k[non_zero_theta[ok_log[num_ok_log - 1]]] * lambda * log(x) / x }  
              
              
              U_1 <- net_stat$Sum_m_k[update_theta] - lambda * (- 2 * extract_u(update_k, update_k) * w_k[update_theta] * (
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
              
              U_2 <- (4 * extract_u(update_k, update_k)^2 * w_k[update_theta] + 8 * extract_u(update_k - 1, update_k)^2 * w_k[minus_1] + 
                        8 * extract_u(update_k + 1, update_k)^2 * w_k[plus_1]) * lambda
              U_3 <- temp4[ok_log[-c(1,2,num_ok_log - 1,num_ok_log)]]
              
            
              
              theta[non_zero_theta[ok_log[1]]] <- tryCatch(uniroot(g_1_new,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                           error = function(e) {return(theta[non_zero_theta[ok_log[1]]])})
              
              theta[non_zero_theta[ok_log[2]]] <- tryCatch(uniroot(g_2,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                           error = function(e) {return(theta[non_zero_theta[ok_log[2]]])})
              #parallelization here
              for (jj in 1:length(update_theta)) {
              #jj <- length(update_theta)
              g <- function(x){U_1[jj]/x - U_2[jj]*log(x)/x - U_3[jj]}
               theta[update_theta[jj]] <- tryCatch(uniroot(g,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                   error = function(e) theta[update_theta[jj]])
              }
               theta[non_zero_theta[ok_log[num_ok_log - 1]]] <- tryCatch(uniroot(g_semiend,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                                        error = function(e) 
                                                                     {return(theta[non_zero_theta[ok_log[num_ok_log - 1]]])})
               theta[non_zero_theta[ok_log[num_ok_log]]] <- tryCatch(uniroot(g_end,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                                    error = function(e) {return(theta[non_zero_theta[ok_log[num_ok_log]]])})
               
               if (length(not_ok_log) > 0) {
                  #print("Update zero")
                  if (lambda <= 0) {
                      theta[non_zero_theta[not_ok_log]] <- net_stat$Sum_m_k[non_zero_theta[not_ok_log]]/temp4[not_ok_log] 
                  }else {
                    g_0  <- function(x) {
                      (net_stat$Sum_m_k[non_zero_theta[not_ok_log]] - 2 * w_k[non_zero_theta[not_ok_log]] * lambda * 
                         (- log(theta[non_zero_theta[not_ok_log]]) - 
                            log(theta[non_zero_theta[ok_log[1]]]) )) / x -   
                        temp4[not_ok_log] - 4 * lambda * w_k[non_zero_theta[not_ok_log]] * log(x) / x}
                    theta[non_zero_theta[not_ok_log]] <- tryCatch(uniroot(g_0,interval = c(0.0000001,1000),tol = .Machine$double.eps)$root,
                                                                  error = function(e) {return(theta[non_zero_theta[not_ok_log]])}) 
                  }
                  
               }
               
               }
              #print(theta)
              #print(cal_reg_A());
            }

         

        } else { # log-linear PA
            .normalized_constant_alpha(normalized_const, alpha,PA_offset,net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
            time_non_zero     <- which(normalized_const != 0)
            non_zero_f        <- which(f != 0) 
            #print(alpha)
            #print(normalized_const[time_non_zero])
            #print(f[non_zero_f])
            
            log_likelihood    <- c(log_likelihood, sum(net_stat$z_j[non_zero_f] * log(f[non_zero_f])) +
                                   alpha * sum(net_stat$Sum_m_k[non_zero_theta] * log(theta[non_zero_theta])) -
                                   sum(net_stat$m_t[time_non_zero] * log(normalized_const[time_non_zero])) + 
                                   ((shape - 1) * (sum(log(f[non_zero_f]))) - rate * sum(f[non_zero_f])) + 
                                   sum(net_stat$offset_m_tk)*log(offset) + (shape - 1) * log(offset) - rate * offset)  
            log_likelihood[length(log_likelihood)] <- log_likelihood[length(log_likelihood)] + net_stat$Sum_m_k[1] * 
                                                   log(PA_offset);
        
            if ((TRUE == debug) && (length(log_likelihood) > 0)){
                print(log_likelihood[length(log_likelihood)])
                if (length(log_likelihood) > 1)
                    if (log_likelihood[length(log_likelihood)] < log_likelihood[length(log_likelihood) - 1])
                        stop("Warning: Log likelihood decreased.")  
            }
            #print("---------")
            break_flag <- FALSE
            if (TRUE == auto_stop)
                if (length(log_likelihood) > 1)
                    tryCatch({if (abs(log_likelihood[length(log_likelihood)] - log_likelihood[length(log_likelihood) - 1]) / 
                               (abs(log_likelihood[length(log_likelihood) - 1]) + 1) < stop_cond)
                    break_flag <- TRUE;},error = function(e) { #print(as.vector(normalized_const));print(f[non_zero_f]);
                  #print(non_zero_f);
                    break_flag <- TRUE;})
           if (break_flag)
               break;   
         #print(theta);
         #flush.console();
        # Remember here theta is always the degree sequence,i.e. we need to calculate the power alpha of this deg. seq.    

          ##################### Update f ######################
        
         if (FALSE == only_PA)  {
             .update_f_alpha(f,non_zero_f,alpha,PA_offset,net_stat$node_degree,theta,net_stat$z_j,normalized_const,net_stat$m_t,shape,rate)
             # update offset
             # if (net_stat$deg_thresh > 0)
             # offset <- .update_offset_alpha(alpha,net_stat$offset_tk, net_stat$offset_m_tk, theta, normalized_const,net_stat$m_t, shape,rate)
         }
         #####################  Update alpha #######################  
         .normalized_constant_alpha(normalized_const, alpha,PA_offset, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
         alpha     <- .update_alpha(non_zero_theta,
                                    normalized_const,f, PA_offset,
                                    theta,net_stat$node_degree,net_stat$m_t,net_stat$Sum_m_k,net_stat$offset_tk,offset) 
         #.normalized_constant_alpha(normalized_const, alpha,PA_offset, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
         #PA_offset <- .update_PA_offset(normalized_const,f,net_stat$node_degree,net_stat$m_t,net_stat$Sum_m_k,
         #                                 net_stat$offset_tk);
         #print(alpha)
         #.normalized_constant_alpha(normalized_const, alpha, PA_offset, net_stat$node_degree,theta,f,net_stat$offset_tk,offset)
          alpha_series <- c(alpha_series,alpha)
      }
}
    theta    <- theta^alpha
    if (only_f == TRUE)
        theta[1] <- PA_offset  
    ########## End of Iteration ########################   

    
    if ((FALSE == break_flag) && (TRUE == auto_stop))
       print(paste0("End by reaching maximum number of iterations (",max_iter,")")); 
    
      if (normalized_f == TRUE) {
          sum_f  <- sum(f)  
          f      <- (length(f) + 1) * f / (sum_f + offset)
          offset <- offset * (length(f) + 1) / (sum_f + offset)   
          if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))){
            .normalized_constant(normalized_const,net_stat$node_degree,theta,f,net_stat$offset_tk,offset) 
          }
          else normalized_const <- as.vector(net_stat$n_tk[,non_zero_theta]%*%theta[non_zero_theta])
      }
       ##### Variance of f ####################
       cov_f   <- rep(0,length(f))
       if (FALSE == only_PA)
           .cal_var_f(cov_f,non_zero_f,net_stat$node_degree,theta,f,net_stat$z_j,
                     normalized_const,net_stat$m_t,shape)
     
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
      hessian_of_regularization_log <- function(theta){
          n      <- length(theta)
          result <- vector()
          if (n - 2 >= 3){
              result <- c(result,2*w_k[non_zero_theta[ok_log[2]]]*((1 - log(theta[1]) / log(center_k[non_zero_theta[ok_log[1]]]) ) /theta[1]^2 + 
                                                                     (2*log(theta[2]) / log(center_k[non_zero_theta[ok_log[2]]]) - 
                                                                        log(theta[3])) / log(center_k[non_zero_theta[ok_log[3]]]) / theta[1]^2 ) )
            
              result <- c(result,2*w_k[non_zero_theta[ok_log[3]]]*((1 - log(theta[2]) / log(center_k[non_zero_theta[ok_log[2]]]) )/theta[2]^2 + 
                                                                   (2*log(theta[3]) / log(center_k[non_zero_theta[ok_log[3]]]) - log(theta[4]) / log(center_k[non_zero_theta[ok_log[4]]])) / theta[2]^2 ) +
                          2*w_k[non_zero_theta[2]]*((2 - 2 * log(theta[2]) / log(center_k[non_zero_theta[ok_log[2]]]) )/theta[2]^2 + ( log(theta[3]) / log(center_k[non_zero_theta[ok_log[3]]]) + 
                                                                                                                                       log(theta[1]) / log(center_k[non_zero_theta[ok_log[1]]])  ) / theta[2]^2 ))
             for (ii in 3:(n-2))
                 result <- c(result,2*w_k[non_zero_theta[ok_log[ii]]]*((2 - 2 * log(theta[ii]) / log(center_k[non_zero_theta[ok_log[ii]]]) )/theta[ii]^2 + 
                                                                     (log(theta[ii + 1]) / log(center_k[non_zero_theta[ok_log[ii + 1]]]) + 
                                                                      log(theta[ii - 1]) / log(center_k[non_zero_theta[ok_log[ii - 1]]]) ) / theta[ii]^2 ) + 
                             2*w_k[non_zero_theta[ok_log[ii + 1]]]*((1 - log(theta[ii]) / log(center_k[non_zero_theta[ok_log[ii]]]) )/theta[ii]^2 + 
                                                                    (2 * log(theta[ii + 1] / log(center_k[non_zero_theta[ok_log[ii + 1]]])) - 
                                                                   log(theta[ii + 2]) / log(center_k[non_zero_theta[ok_log[ii + 2]]]) ) / theta[ii]^2 ) +
                             2*w_k[non_zero_theta[ok_log[ii-1]]]*((1 - log(theta[ii]) / log(center_k[non_zero_theta[ok_log[ii]]]) )/theta[ii]^2 + 
                                                                  (2 * log(theta[ii - 1]) / log(center_k[non_zero_theta[ok_log[ii - 1]]]) - 
                                                                   log(theta[ii - 2]) / log(center_k[non_zero_theta[ok_log[ii - 2]]]) ) / theta[ii]^2 ))
                 ii <- n - 1
                 result <- c(result, 2 * w_k[non_zero_theta[ok_log[ii]]]*((2 - 2 * log(theta[ii]) / log(center_k[non_zero_theta[ok_log[ii]]]) )/theta[ii]^2 + 
                                                                          (log(theta[ii + 1]) / log(center_k[non_zero_theta[ok_log[ii + 1]]]) + 
                                                                           log(theta[ii - 1]) / log(center_k[non_zero_theta[ok_log[ii - 1]]])  ) / theta[ii]^2 ) +
                            2*w_k[non_zero_theta[ok_log[ii-1]]]*((1 - log(theta[ii]) / log(center_k[non_zero_theta[ok_log[ii]]])) / theta[ii]^2 + 
                                                                   (2 * log(theta[ii - 1]) / log(center_k[non_zero_theta[ok_log[ii - 1]]]) - 
                                                                    log(theta[ii - 2]) / log(center_k[non_zero_theta[ok_log[ii - 2]]]) ) / theta[ii]^2 ))
                 ii <- n
                 result <- c(result, 2*w_k[non_zero_theta[ok_log[ii-1]]]*((1 - log(theta[ii]) / log(center_k[non_zero_theta[ok_log[ii]]]) )/theta[ii]^2 + 
                                                                        (2 * log(theta[ii - 1]) / log(center_k[non_zero_theta[ok_log[ii - 1]]]) - 
                                                                           log(theta[ii - 2]) / log(center_k[non_zero_theta[ok_log[ii - 2]]])  ) / theta[ii]^2 )) 
          }
          return(w_k[non_zero_theta[ok_log]]*result)
      }

      # mode_reg_A = 2:
      hessian_of_regularization_mode_2 <- function(theta){
        n      <- length(theta)
        result <- vector()
        if ((n >= 2) && (num_ok_log >= 3)) {
            if (length(not_ok_log) > 0) {  
          result <- c(result, 2 * w_k[non_zero_theta[not_ok_log]] * (1 - log(theta[non_zero_theta[1]]) + log(theta[non_zero_theta[ok_log[1]]])) / theta[non_zero_theta[1]]^2) 
          
          result <- c(result, - (-2 *  w_k[non_zero_theta[not_ok_log]] * log(theta[non_zero_theta[not_ok_log]]) - 2 *  w_k[non_zero_theta[not_ok_log]]  + 
                                  2 *  w_k[non_zero_theta[not_ok_log]] * log(theta[non_zero_theta[ok_log[1]]])  +  
                                  2 *  w_k[non_zero_theta[ok_log[2]]] * u[2,1] * (u[2,3] * log(theta[non_zero_theta[ok_log[3]]]) -  u[2,2] * log(theta[non_zero_theta[ok_log[2]]])) - 
                                  2 *  w_k[non_zero_theta[ok_log[2]]] * u[2,1]^2 + 2 *  w_k[non_zero_theta[ok_log[2]]] * u[2,1]^2 * log(theta[non_zero_theta[ok_log[1]]]) ) / theta[non_zero_theta[ok_log[1]]]^2);
            } else {
              result <- c(result, 0) 
              
              result <- c(result, - (  2 *  w_k[non_zero_theta[not_ok_log]] * log(theta[non_zero_theta[ok_log[1]]])  +  
                                       2 *  w_k[non_zero_theta[ok_log[2]]] * u[2,1] * (u[2,3] * log(theta[non_zero_theta[ok_log[3]]]) -  u[2,2] * log(theta[non_zero_theta[ok_log[2]]])) - 
                                       2 *  w_k[non_zero_theta[ok_log[2]]] * u[2,1]^2 + 2 *  w_k[non_zero_theta[ok_log[2]]] * u[2,1]^2 * log(theta[non_zero_theta[ok_log[1]]]) ) / theta[non_zero_theta[ok_log[1]]]^2); 
           }
          
          result <- c(result, - (-2 *  w_k[non_zero_theta[ok_log[2]]] * (u[2,3] * log(theta[non_zero_theta[ok_log[3]]]) +  u[2,1] * log(theta[non_zero_theta[ok_log[1]]]))  * u[2,2] - 
                                  2 *  w_k[non_zero_theta[ok_log[2]]] * u[2,2]^2 + 2 *  w_k[non_zero_theta[ok_log[2]]] * u[2,2]^2 *  log(theta[non_zero_theta[ok_log[2]]]) + 
                                  2 *  w_k[non_zero_theta[ok_log[3]]] * (u[3,4] * log(theta[non_zero_theta[ok_log[4]]]) -  u[3,3] * log(theta[non_zero_theta[ok_log[3]]]))  * u[3,2] - 
                                  2 *  w_k[non_zero_theta[ok_log[3]]] * u[3,2]^2 +  2 *  w_k[non_zero_theta[ok_log[3]]] * u[3,2]^2 * log(theta[non_zero_theta[ok_log[2]]]) ) / theta[non_zero_theta[ok_log[2]]]^2);
          m <- length(ok_log)
          for (ii in 3:(m-2))
            result <- c(result, - ( 2 *  w_k[non_zero_theta[ok_log[ii - 1]]] * ( - u[ii - 1,ii] + log(theta[non_zero_theta[ok_log[ii]]]) + u[ii - 1,ii - 2] * log(theta[non_zero_theta[ok_log[ii - 2]]]) -
                                                                                   u[ii - 1, ii - 1] * log(theta[non_zero_theta[ok_log[ii - 1]]]))  * u[ii - 1,ii]  + 
                                    2 *  w_k[non_zero_theta[ok_log[ii + 1]]] * ( u[ii + 1,ii + 2] *  log(theta[non_zero_theta[ok_log[ii + 2]]]) -  u[ii + 1,ii] + log(theta[non_zero_theta[ok_log[ii]]]) -
                                                                                  u[ii + 1,ii + 1] * log(theta[non_zero_theta[ok_log[ii + 1]]]))  * u[ii + 1,ii] -    
                                    2 *  w_k[non_zero_theta[ok_log[ii]]] * ( u[ii,ii + 1] *  log(theta[non_zero_theta[ok_log[ii + 1]]]) +  u[ii,ii] - log(theta[non_zero_theta[ok_log[ii]]]) +
                                                                             u[ii,ii - 1] * log(theta[non_zero_theta[ok_log[ii - 1]]]))  * u[ii,ii]  ) / theta[non_zero_theta[ok_log[ii]]]^2)
                        
          ii <- m - 1
          result <- c(result, - ( - 2 *  w_k[non_zero_theta[ok_log[ii]]] * (u[ii,ii + 1] * log(theta[non_zero_theta[ok_log[ii + 1]]]) + u[ii,ii - 1] * log(theta[non_zero_theta[ok_log[ii - 1]]]) +
                                                                                 u[ii, ii] - log(theta[non_zero_theta[ok_log[ii]]]))  * u[ii,ii]  + 
                                    2 *  w_k[non_zero_theta[ok_log[ii - 1]]] * ( - u[ii - 1, ii]  +  log(theta[non_zero_theta[ok_log[ii]]]) + u[ii - 1,ii - 2] *  log(theta[non_zero_theta[ok_log[ii - 2]]]) -  
                                                                                   u[ii - 1,ii - 1] * log(theta[non_zero_theta[ok_log[ii - 1]]]) ) * u[ii - 1,ii]  ) / theta[non_zero_theta[ok_log[ii]]]^2) 
          ii <- m
          result <- c(result, - (2 *  w_k[non_zero_theta[ok_log[ii - 1]]] * (- u[ii - 1, ii] +  log(theta[non_zero_theta[ok_log[ii]]]) + u[ii - 1,ii - 2] * log(theta[non_zero_theta[ok_log[ii - 2]]]) - 
                                                                                  u[ii - 1,ii - 1] * log(theta[non_zero_theta[ok_log[ii - 1]]]) ) * u[ii - 1,ii] / theta[non_zero_theta[ok_log[ii]]]^2 )) 
        }
        return(result);
        
      }

      #interpolation for PA function
      #but only if (FALSE == only_f)
    if ((FALSE == only_f) && (mode_f[1] != "Log_linear")) {  
        if (TRUE == interpolate) {
            theta_nonzero <- which(net_stat$Sum_m_k != 0)
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
                    if (theta_nonzero[i+1] > theta_nonzero[i] + 1) {
                        regress <- lm(c(log(theta[theta_nonzero[i]]),log(theta[theta_nonzero[i+1]]))~ 
                                      c(log(center_k[theta_nonzero[i]]),log(center_k[theta_nonzero[i+1]]))) 
                    for (j in (theta_nonzero[i] + 1):(theta_nonzero[i+1] - 1))
                        theta[j] <- exp(log(center_k[j])*regress$coefficients[2] + 
                                        regress$coefficients[1])           
                    }  
            }
            if (length(theta_nonzero) > 0)
                if (theta_nonzero[length(theta_nonzero)] < length(theta))
                    theta[(theta_nonzero[length(theta_nonzero)] + 1):length(theta)] <- 
                    theta[theta_nonzero[length(theta_nonzero)]]  
        }
    }

    if (only_PA == FALSE) {
        beg     <- which(center_k >= net_stat$deg_thresh & theta != 0)[1]
    } else beg <- which(theta != 0)[1]    
    theta <- theta/theta[beg]
   

    if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))){
        .normalized_constant(normalized_const,net_stat$node_degree,theta,f,net_stat$offset_tk,offset) 
    }
        else normalized_const <- as.vector(net_stat$n_tk[,non_zero_theta]%*%theta[non_zero_theta])
    ############# Calculate the variance of theta #####################
    cov_bin1       <- rep(0,length(theta)) 
    #non_zero       <- which(theta != 0)
    non_zero <- which(net_stat$Sum_m_k > 0)
    non_zero_theta <- non_zero
    time_non_zero     <- which(normalized_const != 0)
    if ((FALSE == only_PA) || ((TRUE == only_PA) && (!is.null(true_f)))) {
      temp4  <- .coeff_var(net_stat$node_degree, f, normalized_const,net_stat$m_t,net_stat$offset_tk, net_stat$start_deg + net_stat$G) 
      temp4 <- temp4[non_zero]
    }
   
    else {
      temp4 <- colSums(net_stat$n_tk[time_non_zero,non_zero_theta]^2 * net_stat$m_t[time_non_zero] / normalized_const[time_non_zero]^2)
    }
    if (0 == mode_reg_A) {
        aa <- 1 / (net_stat$Sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2 + - temp4 + 
                   lambda * hessian_of_regularization(theta[non_zero_theta]))
        bb <-  1 / (net_stat$Sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2 +  
                   lambda * hessian_of_regularization(theta[non_zero_theta])) 
    } else if (1 == mode_reg_A) {
    # mode_reg_A = 1
       
        upper_aa <- net_stat$Sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2 + - temp4 
        upper_bb <- net_stat$Sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2
        
        if (length(ok_log) > 0) {
            upper_aa[ok_log] <- upper_aa[ok_log] + lambda * hessian_of_regularization_log(theta[non_zero_theta[ok_log]])
            upper_bb[ok_log] <- upper_bb[ok_log] + lambda * hessian_of_regularization_log(theta[non_zero_theta[ok_log]])
        }
       
        
        aa <- 1 / (upper_aa)
        bb <-  1 / (upper_bb)   
    } else if (2 == mode_reg_A) {
        upper_aa <- net_stat$Sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2 + - temp4 
        #print(upper_aa)
        #print(non_zero_theta)
        #print(hessian_of_regularization_mode_2(theta))
        temp_aa <- lambda * hessian_of_regularization_mode_2(theta)
        #print(temp_aa)
        if (length(temp_aa) > 0)
            upper_aa <- upper_aa + temp_aa
         
        upper_bb <- net_stat$Sum_m_k[non_zero_theta] /theta[non_zero_theta] ^ 2
        if (length(temp_aa) > 0)
            upper_bb <- upper_bb + lambda * hessian_of_regularization_mode_2(theta)
        aa <- 1 / (upper_aa)
        bb <-  1 / (upper_bb)    
    }
    
     cov_bin1[non_zero_theta] <- ifelse(aa > 10^-10, aa, bb) 
     #print(cov_bin1)
     cov_bin1[cov_bin1 == Inf] <- 0

     cov_bin      <- cov_bin1
     var_log_bin  <- cov_bin / ifelse(theta != 0 , theta ^ 2 , 1)
     upper_bin    <- exp(log(theta) + 2 * sqrt(var_log_bin))
     lower_bin    <- exp(log(theta) - 2 * sqrt(var_log_bin))
     

     non_zero_center <- center2_k > 0 & theta > 0 &  var_log_bin > 0 
     #non_zero_center <- theta > 0 
     if ((FALSE == only_f) && (mode_f[1] != "Log_linear"))
         alpha_center <- lm(log(theta[non_zero_center]) ~ log(center2_k[non_zero_center]))$coefficients[2]
     else 
         alpha_center <- NULL   
     ############# Return theta to A #####################################
      A                                  <- rep(0 , net_stat$deg.max)
      cov                                <- rep(0 , net_stat$deg.max)   
      weight_A                           <- rep(1 , net_stat$deg.max)
      A[1 : (net_stat$start_deg + 1)]        <- theta[1:(net_stat$start_deg + 1)] 
      cov[1 : (net_stat$start_deg + 1)]      <- cov_bin[1:(net_stat$start_deg + 1)]
      weight_A[1 : (net_stat$start_deg + 1)] <- 1
      for (i in 1 : net_stat$G) {
          weight_A[(net_stat$begin_deg[i] : net_stat$end_deg[i]) + 1] <- net_stat$interval_length[i]             
          A[(net_stat$begin_deg[i] : net_stat$end_deg[i]) + 1]        <- theta[net_stat$start_deg + i]
          cov[(net_stat$begin_deg[i] : net_stat$end_deg[i])+1]        <- cov_bin[net_stat$start_deg + i] #* weight_A[net_stat$begin_deg[i] + 1]
      }
      interval     <- 0 : (length(A) - 1)
      
      #return(list(k = interval, A = A))
      non_zero     <- which(A > 10^-20 & cov > 10^-20)
      non_zero     <- non_zero[non_zero >= net_stat$deg_thresh]
      k_non_zero   <- interval[non_zero]
    
      A            <- A[non_zero] 
      #cc           <- exp(mean(log(k_non_zero[k_non_zero > 0])) - mean(log(A[A > 0])))
      cc           <- 1
      A            <- cc * A
      weight_A     <- weight_A[non_zero]
      cov          <- cc ^ 2 * cov[non_zero]  
      ############### fitting A_k = k^alpha ##################
      var_log     <- cov / (A ^ 2)
      sd_log      <- sqrt(var_log)
      log_A       <- log(A)
      log_k       <- log(k_non_zero)
      
      if (((only_f == FALSE) && (mode_f[1] != "Log_linear")) && (length(k_non_zero) > 0))  {
              if (0 == k_non_zero[1]) 
                  linear_fit  <- lm(log_A[-1] ~ log_k[-1] , weights = 1 / (weight_A[-1] * var_log[-1])) 
              else
                  linear_fit  <- lm(log_A ~log_k , weights = 1 / (weight_A * var_log)) 
      }
      else linear_fit <- list(coefficients=c(-1,-1))
      names(linear_fit$coefficients) <- c("offset","Attachment exponent")
      upper_A         <- exp(log(A) + 2 * sd_log)
      lower_A         <- exp(log(A) - 2 * sd_log)
      
     
      ci <- confint(linear_fit,"Attachment exponent")
      
      ##### Variance of f ####################
    
      #if (rate <= 1 & shape <= 1)
      #   f[-non_zero_f] <- 0

      f_new                                    <- rep(offset,net_stat$N)
      names(f_new)                             <- net_stat$node_id
      f_new[as.character(net_stat$f_position)]     <- f
      cov_f_new                                <- rep(0,net_stat$N)
      names(cov_f_new)                         <- net_stat$node_id
      cov_f_new[as.character(net_stat$f_position)] <- abs(cov_f)
      non_zero_f                               <- f_new > 10^-20 & cov_f_new > 10^-20
    
      upper_f                                  <- rep(0,net_stat$N)
      upper_f[non_zero_f]                      <- exp(log(f_new[non_zero_f]) + 2 * sqrt(cov_f_new[non_zero_f] / f_new[non_zero_f] ^ 2))

      lower_f                                  <- rep(0,net_stat$N)
      lower_f[non_zero_f]                      <- exp(log(f_new[non_zero_f]) - 2 * sqrt(cov_f_new[non_zero_f] / f_new[non_zero_f] ^ 2))

      if (mode_f[1] != "Log_linear")
         alpha = linear_fit$coefficients[2]
      else names(alpha) <- "Estimated attachment exponent"
      
      result <- list(# estimated PA function and its variances, confidence interval 
                     k       = k_non_zero ,  A             = A          , var_A       = cov       , var_logA = var_log,
                     upper_A = upper_A    ,  lower_A       = lower_A    , weight_of_A = weight_A  , center_k = center_k,  
                     theta   = theta      ,  upper_bin     = upper_bin  , lower_bin   = lower_bin , var_bin  = cov_bin,
                     # estimated attachment exponent alpha, and the log-linear fit 
                     alpha   = alpha      ,  loglinear_fit = linear_fit , ci          = ci        ,
                     
                     # if mode_f = "Log_linear", the  attachment exponent alpha's over iterations
                     alpha_series = ifelse(rep(mode_f[1] == "Log_linear",length(alpha_series)),alpha_series,-1),
                     
                     # estimated node fitnesses and their variances, confidence intervals
                     f        = f_new     ,  var_f         = cov_f_new, 
                     upper_f  = upper_f   ,  lower_f       = lower_f, # confidence intervals
                     
                     # values of the objective function over iterations
                     objective_value = log_likelihood,
                     
                     # other parameters specified
                     mode_f = mode_f[1] , true_A = true_A , true_f = true_f, PA_offset = PA_offset, candidate_accept = candidate_accept,
                     only_PA = only_PA, only_f = only_f, lambda = lambda, shape = shape, rate = rate, normalized_f = normalized_f, 
                     deg_threshold = net_stat$deg_thresh, stop_cond = stop_cond, auto_lambda = auto_lambda, ratio = ratio, 
                     G = net_stat$G,shape = shape, rate = rate, offset = offset)
      class(result) <- "PAFit_result"
      return(result)
}
