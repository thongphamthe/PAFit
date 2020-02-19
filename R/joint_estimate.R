joint_estimate <- function(net_object                                 , 
                           net_stat       = get_statistics(net_object), 
                           p              = 0.75                      ,
                           stop_cond      = 10^-8                     ,
                           mode_reg_A     = 0                         ,
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
  
  # first finding the optima r and s by cross validation
  
  deg_thresh     <- net_stat$deg_thresh
  cv_deg_thresh  <- c(1);
  
  #net_type       <- net_stat$net_type
  #raw_net_object <- as.PAFit_net(graph = raw_net, type = net_type)
  
  data_cv   <- .CreateDataCV(net_object, g = net_stat$g, deg_thresh = deg_thresh,
                             p = p)
  
  cv_result <- .performCV_old_new(data_cv, stop_cond = stop_cond, mode_reg_A = mode_reg_A, 
                                 print_out = FALSE, cv_deg_thresh = cv_deg_thresh, 
                                 normal_start_f = TRUE, weight_f = 0)
  
  
  
  s_optimal      <- cv_result$s_optimal
  r_optimal      <- cv_result$r_optimal
  lambda_optimal <- cv_result$lambda_optimal
  
  alpha_optimal  <- cv_result$alpha_optimal
  #print(s_optimal)
  #print(r_optimal)
  
  
  f_vector        <- rep(1,length(net_stat$f_position))
  names(f_vector) <- net_stat$f_position
  
  name_vec <- as.character(net_stat$f_position)
  name_cv  <- names(cv_result$estimated_fitness)
  for (ii in 1:length(name_vec)) {
    jj <- which(name_cv == name_vec[ii])  
    if (length(jj) > 0)
      f_vector[name_vec[ii]] <- cv_result$estimated_fitness[jj]
  }
  
  #print(f_vector)
  
  # find a rough estimate of node fitnesses and attachment function based on the model Ak = k^alpha
  # use the estimated alpha and node fitnesses of the learning data as warm-start position
  
  #result_temp     <- PAFit(net_stat,
  #                         mode_f       = "Log_linear",
  #                         s            = cv_result$s_optimal,
  #                         start_f      = f_vector,
  #                         stop_cond    = stop_cond * 100 # loose convergence condition
  #)
  # feed the estimated attachment function and node fitnesses for a warm-start rerun with nonparametric attachment function
  
  result  <- PAFit(net_stat,
                   r            = cv_result$r_optimal,
                   s            = cv_result$s_optimal,
                   start_f      = f_vector           ,
                   alpha_start  = alpha_optimal      ,
                   stop_cond    = stop_cond          ,
                   weight_power = 0                  ,
                   mode_reg_A   = mode_reg_A)
  
  
  combined_result <- list(cv_data = data_cv, cv_result = cv_result, 
                          estimate_result = result) 
  class(combined_result) <- "Full_PAFit_result"
  return(combined_result)
}