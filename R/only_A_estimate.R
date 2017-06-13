only_A_estimate <- function(net_object         , 
                            net_stat           , 
                            stop_cond  = 10^-8 , 
                            mode_reg_A = 0     ,
                            ...) {
  
  if (class(net_object) != "PAFit_net")
    stop("net_object should be of PAFit_net class.")
  
  if (class(net_stat) != "PAFit_data")
    stop("Please input a proper net summary of class PAFit_data");
  
  net_type       <- net_stat$net_type

  data_cv        <- .CreateDataCV_onlyA(net_object, deg_thresh = net_stat$deg_thresh)
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
  combined_result        <- list(cv_data = data_cv, cv_result = cv_result, estimate_result = result) 
  class(combined_result) <- "Full_PAFit_result"
  return(combined_result)
}