only_F_estimate <- function(net_object, 
                            net_stat  = get_statistics(net_object), 
                            p         = 0.75      ,
                            stop_cond = 10^-8     , 
                            model_A   = "Linear"  ,
                            ...) {
  
  if (!is(net_object,"PAFit_net"))
    stop("net_object should be of PAFit_net class.")
  
  if (!is(net_stat,"PAFit_data"))
    stop("Please input a proper net summary of class PAFit_data");
  
  mode_f <- "Log_linear"
  if (model_A == "Linear")
      alpha_start <- 1    
  else if (model_A == "Constant")
      alpha_start   <- 0  
  
  net_type     <- net_stat$net_type
  
  data_cv        <- .CreateDataCV(net_object, deg_thresh = net_stat$deg_thresh, p = p)
  cv_result      <- .OnlyF_CV(data_cv, 
                            stop_cond   = stop_cond, 
                            alpha       = alpha_start,
                           ...)
  f_vector        <- rep(1,length(net_stat$f_position))
  #print(paste0("length f_vector before: ",length(f_vector)))
  names(f_vector) <- net_stat$f_position
  #print(length(cv_result$estimated_fitness))
  f_vector[names(cv_result$estimated_fitness)] <- cv_result$estimated_fitness
  #print(paste0("length f_vector after: ",length(f_vector)))
  
  result  <- PAFit(net_stat, 
                   only_f      = TRUE                 , 
                   mode_f      = mode_f               , 
                   start_f     = f_vector             ,
                   s           = cv_result$s_optimal  , 
                   alpha_start = alpha_start          ,
                   stop_cond   = stop_cond            , 
                   ...) 
  #stop("reach here")
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
  
  class(combined_result) <- "Full_PAFit_result"
  return(combined_result)
}