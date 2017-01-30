Jeong <- function(raw_net, net_stat,T_0,T_1, interpolate = FALSE) {
    #check the class 
    if (is.null(raw_net))
        stop("Please input a 3-column matrix (see references for more information)");
    if (is.null(net_stat))
        stop("Please input a proper net summary (see the class pafit_summary in references for more information)");
    if (is.null(T_0))
       stop("Please input a positive integer for the starting time T_0")
    if (is.null(T_1))
        stop("Please input a positive integer for the end time T_1")
    if ((0 >= T_0) || (0 >= T_1))
        stop("The time must be positive integers")  
   
    unique_time   <- sort(unique(raw_net[,3]))
    T_0           <- unique_time[T_0]
    T_1           <- unique_time[T_1]
    if (T_1 > length(unique_time))
        stop("The ending time T_1 is too large. Its maximum value is the number of distinct time-stamps of the data.")  
    in_node       <- raw_net[,2]
    out_node      <- raw_net[,1]
    time          <- raw_net[,3]
    N_0           <- sort(union(in_node[time <= T_0],out_node[time <= T_0]))
    degree        <- rep(0,length(N_0))
    names(degree) <- N_0
    count_deg     <- table(in_node[time <= T_0])
    degree[labels(count_deg)[[1]]] <- net_stat$bin_vector[count_deg]
  
    to_node   <- in_node[(time > T_0) & (time <= T_1)]
    from_node    <- out_node[time > T_0 & time <= T_1]
  N_1       <- setdiff(union(to_node,from_node),N_0)
  
  #if (net_stat$Sum_m_k[length(net_stat$Sum_m_k)] != 0)
      count            <- rep(0,net_stat$G + net_stat$start_deg)
  #else 
  #    count            <- rep(0,net_stat$G + net_stat$start_deg - 1)
  
  normalized_const <- rep(0,length(count))
  
  temp      <- table(degree)
  for (i in 1:length(degree))
    normalized_const[as.integer(labels(temp)[[1]]) + 1] <- temp
  
  really <- which((from_node %in% N_1 & to_node %in% N_0) == TRUE)
  temp   <- table(degree[as.character(to_node[really])])
  count[as.integer(labels(temp)[[1]]) + 1] <- temp
  
  count <- count/ifelse(normalized_const != 0,normalized_const,1) 
  
  theta <- count
  

  
  center_k  <- rep(0, length(theta))
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
  
  count <- theta
  if (net_stat$Binning == TRUE) {
      A                             <- rep(0,net_stat$deg.max)
      A[1:(net_stat$start_deg + 1)] <- theta[1:(net_stat$start_deg + 1)] 
      for (i in 1:net_stat$G) {
          A[(net_stat$begin_deg[i]:net_stat$end_deg[i]) + 1]        <- theta[net_stat$start_deg + i]
      }
  } else 
      A <- theta
  
  if (net_stat$Sum_m_k[length(net_stat$Sum_m_k)] == 0) 
      A <- A[-length(A)]
  
  A[which(A == "NaN")] <- 0
  k <- 0:(length(A) - 1)
  

  
  result        <- list(A = A , k = k,center_k = center_k, theta = count, alpha = alpha, loglinear_fit = linear_fit) 
  class(result) <- "PA_result"
  return(result)
}