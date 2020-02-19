Jeong <- function(net_object                              , 
                  net_stat    = get_statistics(net_object), 
                  T_0_start   = 0                         ,
                  T_0_end     = round(net_stat$T * 0.75) ,
                  T_1_start   = T_0_end + 1               ,
                  T_1_end     = net_stat$T                ,
                  interpolate = FALSE) {
    #check the class 
    if (!is(net_object,"PAFit_net"))
      stop("net_object should be of PAFit_net class.")
    
    raw_net <- net_object$graph
  
    if (!is(net_stat,"PAFit_data"))
        stop("Please input a proper net summary of class PAFit_data");
    
    if (T_0_end < T_0_start) 
        stop(cat("T_0_end must be at least T_0_start"))    
  
    if (T_1_start <= T_0_end) 
        stop(cat("T_1_start must be greater than T_0_end"))    
  
    if (T_1_end < T_1_start) 
        stop(cat("T_1_end must be at least T_1_start"))
  
    unique_time   <- sort(unique(raw_net[,3]))
    if (T_1_end > length(unique_time))
        stop("The ending interval T_1 is too large. Its maximum value is the number of distinct time-stamps of the data.")  
    
    T_0_interval  <- unique_time[unique_time >= T_0_start & unique_time <= T_0_end]
    T_1_interval  <- unique_time[unique_time >= T_1_start & unique_time <= T_1_end]
    
    in_node       <- raw_net[,2]
    out_node      <- raw_net[,1]
    time          <- raw_net[,3]
    
    N_0           <- sort(union(in_node[time %in% T_0_interval],out_node[time %in% T_0_interval]))
    degree        <- rep(0,length(N_0))
    names(degree) <- N_0
    
    # freeza the degrees of N_0 nodes at the end of interval T_0_interval
    count_deg     <- table(in_node[time %in% T_0_interval])
    degree[labels(count_deg)[[1]]] <- net_stat$bin_vector[count_deg]
  
    to_node       <- in_node[time %in% T_1_interval]
    from_node     <- out_node[time %in% T_1_interval]
    
   #if (net_stat$sum_m_k[length(net_stat$sum_m_k)] != 0)
    count            <- rep(0,net_stat$g + net_stat$start_deg)
  #else 
  #    count            <- rep(0,net_stat$g + net_stat$start_deg - 1)
  
    normalized_const <- rep(0,length(count))
  
    temp      <- table(degree)
    for (i in 1:length(degree))
        normalized_const[as.integer(labels(temp)[[1]]) + 1] <- temp
  
    really <- which((to_node %in% N_0) == TRUE)
    temp   <- table(degree[as.character(to_node[really])])
    count[as.integer(labels(temp)[[1]]) + 1] <- temp
  
    count <- count/ifelse(normalized_const != 0,normalized_const,1) 
  
  theta <- count
  
  
  
  # center_k  <- rep(0, length(theta))
  # if (net_stat$start_deg > 0)
  #   center_k[1:net_stat$start_deg]  <- 0:(net_stat$start_deg - 1)
  # for (i in 1:net_stat$g) {
  #   if (net_stat$begin_deg[i] != 0) {
  #     #              center_k[i]  <- round((net_stat$begin_deg[i] + net_stat$end_deg[i])/2)  
  #     center_k[i] <- round(net_stat$begin_deg[i]*sqrt((net_stat$begin_deg[i] + net_stat$interval_length[i] - 1)/ net_stat$begin_deg[i]))
  #   } else
  #     center_k[i] <- net_stat$end_deg[i]  
  #   #          center_k[i]  <- round((net_stat$begin_deg[i] + net_stat$end_deg[i])/2)   
  #   #        
  # }
  center_k           <- net_stat$center_k
  center_k_degthresh <- which(center_k == net_stat$deg_thresh)[1]
  if (length(center_k_degthresh) > 0) {
      if (!is.na(theta[center_k_degthresh]))    
          if (theta[center_k_degthresh] !=0 ) {
              theta <- theta / theta[center_k_degthresh]  
      }
  }
  
  #estimate the attachment exponent alpha
  non_zero     <- which(center_k > 0 & theta > 0)
  log_A        <- log10(theta[non_zero])
  log_k        <- log10(center_k[non_zero])
  linear_fit   <- lm(log_A ~ log_k)
  alpha        <- linear_fit$coefficients[2]
  names(alpha) <- "Estimated attachment exponent"
  names(linear_fit$coefficients) <- c("Constant","Attachment exponent")
  res          <- df.residual(linear_fit)
  if (res > 0)
      ci       <- confint(linear_fit,"Attachment exponent")
  else ci      <- "N" 
  

  if (TRUE == interpolate) {
      theta_nonzero <- which(theta != 0)
      if (length(theta_nonzero) > 0)
          if (theta_nonzero[1] > 1) {
              theta[1:(theta_nonzero[1] - 1)] <- theta[theta_nonzero[1]]  
          }
      if (length(theta_nonzero) > 1) {
          for (i in 1:(length(theta_nonzero) - 1))
              if (theta_nonzero[i+1] > theta_nonzero[i] + 1) {
                  regress_flag <- 0  
                  if (center_k[theta_nonzero[i]] > 0 && center_k[theta_nonzero[i+1]] > 0 &&
                      theta[theta_nonzero[i]] > 0 && theta[theta_nonzero[i+1]] > 0){
                          regress_flag <- 1
                          regress <- lm(c(log10(theta[theta_nonzero[i]]),log10(theta[theta_nonzero[i+1]]))~ 
                                      c(log10(center_k[theta_nonzero[i]]),log10(center_k[theta_nonzero[i+1]])))
                  } else if (center_k[theta_nonzero[i]] > 0 && center_k[theta_nonzero[i+1]] > 0 &&
                                theta[theta_nonzero[i]] > 0 && theta[theta_nonzero[i+1]] > 0) {
                             regress_flag <- 1
                             regress <- lm(c(log10(theta[theta_nonzero[i]]) , log10(theta[theta_nonzero[i+1]]))~ 
                                        c(log10(center_k[theta_nonzero[i]] + 1) , log10(center_k[theta_nonzero[i+1]] + 1)))
                  }
                 if (1 == regress_flag)
                     for (j in (theta_nonzero[i] + 1) : (theta_nonzero[i+1]-1))
                         theta[j] <- exp(log10(center_k[j]) * regress$coefficients[2] + regress$coefficients[1])           
        }  
    }
    if (length(theta_nonzero) > 0)
      if (theta_nonzero[length(theta_nonzero)] < length(theta))
        theta[(theta_nonzero[length(theta_nonzero)] + 1):length(theta)] <- 
          theta[theta_nonzero[length(theta_nonzero)]] 
  }
  
  count <- theta
  if (net_stat$binning == TRUE) {
      A                             <- rep(0,net_stat$deg_max)
      A[1:(net_stat$start_deg + 1)] <- theta[1:(net_stat$start_deg + 1)] 
      for (i in 1:net_stat$g) {
          A[(net_stat$begin_deg[i]:net_stat$end_deg[i]) + 1]        <- theta[net_stat$start_deg + i]
      }
  } else 
      A <- theta
  
  if (net_stat$sum_m_k[length(net_stat$sum_m_k)] == 0) 
      A <- A[-length(A)]
  
  A[which(A == "NaN")] <- 0
  k <- 0:(length(A) - 1)
  

  result        <- list(A      = A          , k     = k     , center_k      = center_k , 
                        theta  = count      , alpha = alpha , loglinear_fit = linear_fit,
                        g      = net_stat$g ,
                        ci     = ci) 
  class(result) <- "PA_result"
  return(result)
}