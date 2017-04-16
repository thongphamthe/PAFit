only_F_estimate <- function(raw_net, 
                           net_stat, 
                           stop_cond = 10^-9     , 
                           model_A   = "Linear"  ,
                           ...) {
  mode_f <- "Log_linear"
  if (model_A == "Linear")
      alpha_start <- 1    
  else if (model_A == "Constant")
      alpha_start   <- 0  
  
  net_type     <- net_stat$net_type
  
  data_cv      <- .CreateDataCV(net = raw_net, deg_thresh = 0, net_type = net_type)
  cv_result    <- .OnlyF_CV(data_cv, 
                            stop_cond   = stop_cond, 
                            alpha       = alpha_start,
                           ...)
  f_vector        <- rep(1,length(net_stat$f_position))
  names(f_vector) <- net_stat$f_position
  f_vector[names(cv_result$estimated_fitness)] <- cv_result$estimated_fitness
  
  result  <- PAFit(net_stat, 
                   only_f      = TRUE                 , 
                   mode_f      = mode_f               , 
                   start_f     = f_vector             ,
                   s           = cv_result$s_optimal  , 
                   alpha_start = alpha_start          ,
                   stop_cond   = stop_cond            , 
                   ...) 
  combined_result        <- list(cv_data = data_cv, cv_result = cv_result, estimate_result = result) 
  class(combined_result) <- "Full_PAFit_result"
  return(combined_result)
}