JointEstimate <- function(raw_net                  , 
                          net_stat                 , 
                          stop.cond      = 10^-8   ,
                          mode_reg_A     = 0       ,
                          print.out      = FALSE   ,
                          cv_deg_thresh  = c(1,10) ,
                          p              = 0.75    ,
                          ...) {
  
  # first finding the optima r and s by cross validation
  #print("start")
  
  normal_start_f <- 1
  weight_f       <- 0
  deg_thresh     <- net_stat$deg_thresh
  
  #cv_deg_thresh <- c(1,5);
  net_type  <-  net_stat$net_type
  data_cv   <- .CreateDataCV(raw_net, G = net_stat$G, deg_thresh = deg_thresh,
                             p = p, net_type = net_type)
  cv_result <- .performCV_old(data_cv, stop_cond = stop.cond, mode_reg_A = mode_reg_A, 
                              print_out = print.out, cv_deg_thresh = cv_deg_thresh, 
                              normal_start_f = normal_start_f, weight_f = weight_f)
  
  
  s_optimal      <- cv_result$s_optimal
  r_optimal      <- cv_result$r_optimal
  lambda_optimal <- cv_result$lambda_optimal
  
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
                           mode_f       = "Log_linear",
                           s            = cv_result$s_optimal,
                           start_f      = f_vector,
                           alpha_start  = cv_result$alpha,
                           #auto_stop   = FALSE,
                           #iter        = 20,
                           stop_cond    = stop.cond * 100, # loose convergence condition
                           ...)
  # feed the estimated attachment function and node fitnesses for a warm-start rerun with nonparametric attachment function
  #print("Reached final here")
  result  <- PAFit(net_stat,
                   r            = cv_result$r_optimal,
                   #lambda      = lambda_optimal     ,
                   #auto_lambda = FALSE              ,
                   s            = cv_result$s_optimal,
                   start_f      = result_temp$f[as.character(net_stat$f_position)],
                   alpha_start  = result_temp$alpha  ,
                   stop_cond    = stop.cond          ,
                   weight_power = weight_f           ,
                   mode_reg_A   = mode_reg_A         ,
                   ...)
  result_pow_0 <-   result  <- PAFit(net_stat,
                                        r            = cv_result$r_optimal,
                                        #lambda      = lambda_optimal     ,
                                        #auto_lambda = FALSE              ,
                                        s            = cv_result$s_optimal,
                                        start_f      = result_temp$f[as.character(net_stat$f_position)],
                                        alpha_start  = result_temp$alpha  ,
                                        stop_cond    = stop.cond          ,
                                        weight_power = 0                  ,
                                        mode_reg_A   = mode_reg_A         ,
                                        ...) 
  
  combined_result        <- list(cv_data = data_cv, cv_result = cv_result, 
                                 estimate_result = result, estimate_result_pow0 = result_pow_0) 
  class(combined_result) <- "Full_PAFit_result"
  return(combined_result)
}