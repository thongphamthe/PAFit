Newman_corrected <- function(net_stat,start = 1,interpolate = TRUE){
  deg.max          <- net_stat$deg.max
  T_time           <- dim(net_stat$n_tk)[1]      
  k                <- 0:deg.max
  temp             <- net_stat$m_tk[start:T_time,]/net_stat$n_tk[start:T_time,]
  temp[which(temp == "NaN")] <- 0
  temp[which(temp == "Inf")] <- 0
  num_of_node_t    <- rowSums(net_stat$n_tk[start:T_time,])  
  normalized_const <- apply(net_stat$n_tk[start:T_time,],2,function(x) sum(x != 0))
  theta            <- colSums(temp*num_of_node_t)
  theta            <- theta/normalized_const
  #while (0 == theta[length(theta)])
  #    theta <- theta[-length(theta)]
  
  theta[which(theta == "NaN")] <- 0
  theta[which(theta == "Inf")] <- 0 
  center_k  <- rep(0, length(theta))
  theta     <- theta/theta[1]
  if (net_stat$start_deg > 0)
    center_k[1:net_stat$start_deg]  <- 0:(net_stat$start_deg - 1)
  for (i in 1:net_stat$G) {
    if (net_stat$begin_deg[i] != 0) {
      #              center_k[i]  <- round((net_stat$begin_deg[i] + net_stat$end_deg[i])/2)  
      center_k[i] <- round(net_stat$begin_deg[i]*sqrt((net_stat$begin_deg[i] + net_stat$interval_length[i] - 1)/ net_stat$begin_deg[i]))
    } else
      center_k[i] <- net_stat$end_deg[i]  
    #          center_k[i]  <- round((net_stat$begin_deg[i] + net_stat$end_deg[i])/2)   
    #        
  }
  
  #estimate the attachment exponent alpha
  non_zero      <- which(center_k > 0 & theta > 0)
  log_A        <- log(theta[non_zero])
  log_k        <- log(center_k[non_zero])
  
  linear_fit   <- lm(log_A ~ log_k)
  alpha        <- linear_fit$coefficients[2]
  names(alpha) <- "Estimated attachment exponent" 
  
  
  #interpolation
  if (TRUE == interpolate) {
    theta_nonzero <- which(theta != 0)
    if (length(theta_nonzero) > 0)
      if (theta_nonzero[1] > 1) {
        theta[1:(theta_nonzero[1] - 1)] <- theta[theta_nonzero[1]]  
      }
    if (length(theta_nonzero) > 1) {
      for (i in 1:(length(theta_nonzero) - 1))
        if (theta_nonzero[i+1] > theta_nonzero[i] + 1) {
          if (center_k[theta_nonzero[i]] != 0) {
            regress <- lm(c(log(theta[theta_nonzero[i]]),log(theta[theta_nonzero[i+1]]))~ 
                            c(log(center_k[theta_nonzero[i]]),log(center_k[theta_nonzero[i+1]])))
          } else
            regress <- lm(c(log(theta[theta_nonzero[i]]),log(theta[theta_nonzero[i+1]]))~ 
                            c(log(center_k[theta_nonzero[i]] + 1),
                              log(center_k[theta_nonzero[i+1]] + 1)))
          
          for (j in (theta_nonzero[i] + 1):(theta_nonzero[i+1]-1))
            theta[j] <- exp(log(center_k[j]) * regress$coefficients[2] + regress$coefficients[1])           
        }  
    }
    if (length(theta_nonzero) > 0)
      if (theta_nonzero[length(theta_nonzero)] < length(theta))
        theta[(theta_nonzero[length(theta_nonzero)] + 1):length(theta)] <- 
          theta[theta_nonzero[length(theta_nonzero)]] 
  }
  
  
  if (net_stat$Binning == TRUE) {
    
    A                                <- rep(0,net_stat$deg.max)
    #weight_A                         <- rep(0,net_stat$deg.max)  
    #weight_A[1:(net_stat$start_deg + 1)] <- 1
    A[1:(net_stat$start_deg + 1)]       <- theta[1:(net_stat$start_deg + 1)] 
    for (i in 1:net_stat$G) {
      #weight_A[(net_stat$begin_deg[i]:net_stat$end_deg[i]) + 1] <- net_stat$theta_w[net_stat$start_deg + i]
      A[(net_stat$begin_deg[i]:net_stat$end_deg[i]) + 1]        <- theta[net_stat$start_deg + i]
    }
    #A <- A/weight_A  
  }   
  else A <- theta
  if (net_stat$Sum_m_k[length(net_stat$Sum_m_k)] == 0) 
      A <- A[-length(A)]
  
  A[which(A == "NaN")] <- 0
  k                    <- 0:(length(A) - 1)

  
  result        <- list(A = A , k = k,center_k = center_k, theta = theta, alpha = alpha, loglinear_fit = linear_fit) 
  class(result) <- "PA_result"
  
  return(result)
} 