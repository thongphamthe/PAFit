only_A_estimate <- function(net_object                             , 
                            net_stat   = get_statistics(net_object), 
                            p          = 0.75                      , 
                            stop_cond  = 10^-8                     , 
                            mode_reg_A = 0                         ,
                            MLE        = FALSE                     ,
                            ...) {
  
  if (!is(net_object,"PAFit_net"))
    stop("net_object should be of PAFit_net class.")
  
  if (!is(net_stat,"PAFit_data"))
    stop("Please input a proper net summary of class PAFit_data");
  
  # quick check 
  non_zero_theta     <- which(net_stat$sum_m_k > 0)
  num_nonzero        <- length(non_zero_theta)
  if (num_nonzero == 1) {
    # only one non-zero bin
    stop(paste0("Error: There is only one bin that has a non-zero number of new edges (bin ",which(net_stat$sum_m_k > 0),"). To estimate the PA function, we need at least two bins with non-zero number of new edges."))  
  }
  if (num_nonzero == 0) {
    # no non-zero bin
    stop(paste0("Error: There is no bin that has a non-zero number of new edges. To estimate the PA function, we need at least two bins with non-zero number of new edges."))  
  }
  
  if (MLE == FALSE) {
      net_type       <- net_stat$net_type

      data_cv        <- .CreateDataCV_onlyA(net_object, deg_thresh = net_stat$deg_thresh, p = p)
      cv_result      <- .OnlyA_CV(data_cv, stop_cond = stop_cond * 10, mode_reg_A = mode_reg_A,...)
      # find a rough estimate of the attachment function based on the model Ak = k^alpha
      result_temp  <- PAFit(net_stat, 
                            mode_f      = "Log_linear"            , 
                            only_PA     = TRUE                    ,
                            alpha_start = cv_result$alpha_optimal , 
                            stop_cond   = stop_cond * 10          , # loose convergence condition 
                            ...) 
      # feed the estimated attachment function for a warm-start re-run with nonparametric attachment function
      result  <- PAFit(net_stat, 
                   r           = cv_result$r_optimal     , 
                   mode_f      = "Linear_PA"             ,  
                   only_PA     = TRUE                    ,
                   #alpha_start = result_temp$alpha   , 
                   alpha_start = cv_result$alpha_optimal , 
                   stop_cond   = stop_cond               ,
                   mode_reg_A  = mode_reg_A              ,  
                   ...) 
  } else {
        result  <- PAFit(net_stat, 
                         r           =  0, 
                         mode_f      = "Linear_PA"             ,  
                         only_PA     = TRUE                    ,
                         stop_cond   = stop_cond               ,
                         ...)
    data_cv   <- NULL
    cv_result <- NULL
  }
  PA  <- result$A
  fit <- result$f
  small_t <- dim(net_stat$node_degree)[1]
  contrib_PA_array  <- rep(0,small_t)
  contrib_fit_array <- rep(0,small_t) 
  name_node <- colnames(net_stat$node_degree)
  
  for (i in 1:small_t) {
    presence <- net_stat$node_degree[i,] > 0
    
    # sampling the node based on the product of PA and fitness
    pa_value                  <- PA[net_stat$node_degree[i,presence] + 1]
    pa_value[is.na(pa_value)] <- PA[length(PA)]
    #print(length(pa_value))
    
    fitness_value <- fit[name_node[presence]]
    #print(length(fitness_value))
    sampling_prob <- fitness_value * pa_value / sum(fitness_value * pa_value)
    
    mean_log_PA   <- mean(sampling_prob * log(pa_value), na.rm = TRUE)
    var_log_PA    <- mean(sampling_prob* (log(pa_value) - mean_log_PA)^2, na.rm = TRUE)
    
    mean_log_fit   <- mean(sampling_prob * log(fitness_value) , na.rm = TRUE)
    var_log_fit    <- mean(sampling_prob* (log(fitness_value) - mean_log_fit)^2 , na.rm = TRUE)
    
    contrib_PA    <- var_log_PA
    contrib_fit   <- var_log_fit
    contrib_PA_array[i]  <- contrib_PA
    contrib_fit_array[i] <- contrib_fit
  }
  mean_PA_contrib  <- sqrt(mean(contrib_PA_array, na.rm = TRUE))
  mean_fit_contrib <- sqrt(mean(contrib_fit_array, na.rm = TRUE))
  contribution <- list(PA_contribution = sqrt(contrib_PA_array),
                       fit_contribution = sqrt(contrib_fit_array),
                       mean_PA_contrib = mean_PA_contrib,
                       mean_fit_contrib = mean_fit_contrib)
  
  combined_result <- list(cv_data = data_cv, cv_result = cv_result, 
                          estimate_result = result, contribution = contribution) 
  
  class(combined_result) <- "Full_PAFit_result"
  return(combined_result)
}