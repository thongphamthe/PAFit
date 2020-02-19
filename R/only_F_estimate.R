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
  combined_result        <- list(cv_data = data_cv, cv_result = cv_result, estimate_result = result) 
  class(combined_result) <- "Full_PAFit_result"
  return(combined_result)
}