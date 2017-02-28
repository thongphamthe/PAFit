JointEstimate <- function(raw_net, 
                          net_stat, 
                          stop.cond   = 10^-8 ,
                          mode_reg_A  = 0     , 
                          print.out   = FALSE ,
                          ...) {
  
  # first finding the optima r and s by cross validation
  #print("start")
  data_cv   <- .CreateDataCV(raw_net, deg_thresh = 0)
  cv_result <- .performCV(data_cv, stop_cond = stop.cond, mode_reg_A = mode_reg_A, print_out = print.out)
  
 
  s_optimal <- cv_result$s_optimal
  r_optimal <- cv_result$r_optimal
 
  f_vector        <- rep(1,length(net_stat$f_position))
  names(f_vector) <- net_stat$f_position
  #print(length(cv_result$estimated_fitness))
  name_vec <- as.character(net_stat$f_position)
  name_cv  <- names(cv_result$estimated_fitness)
  for (ii in 1:length(name_vec)) {
      jj <- which(name_cv == name_vec[ii])  
      if (length(jj) > 0)
          f_vector[name_vec[ii]] <- cv_result$estimated_fitness[jj]
  }
  
  
  # find a rough estimate of node fitnesses and attachment function based on the model Ak = k^alpha
  # use the estimated alpha and node fitnesses of the learning data as warm-start position
  #print("Reached here")
  #print(f_vector)
  #print(paste0("start f vector: ",length(f_vector)))
  #print(paste0("length f_position (outside): ", length(net_stat$f_position)))
  result_temp     <- PAFit(net_stat,
                           mode_f      = "Log_linear",
                           s           = cv_result$s_optimal,
                           start_f     = f_vector,
                           alpha_start = cv_result$alpha,
                           #auto_stop   = FALSE,
                           #iter        = 20,
                           stop_cond   = stop.cond * 100, # loose convergence condition
                           ...)
  # feed the estimated attachment function and node fitnesses for a warm-start rerun with nonparametric attachment function
  #print("Reached final here")
  result  <- PAFit(net_stat,
                   r           = cv_result$r_optimal,
                   s           = cv_result$s_optimal,
                   start_f     = result_temp$f[as.character(net_stat$f_position)],
                   alpha_start = result_temp$alpha,
                   stop_cond   = stop.cond,
                   mode_reg_A  = mode_reg_A, ...)

  return(list(cv_data = data_cv, cv_result = cv_result, estimate_result = result))
}