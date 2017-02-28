OnlyA_Estimate <- function(raw_net             , 
                           net_stat            , 
                           stop.cond  = 10^-8  , 
                           mode_reg_A = 0      ,
                           ...) {
  
  
  data_cv      <- .CreateDataCV_onlyA(net = raw_net, deg_thresh = 0)
  cv_result    <- .OnlyA_CV(data_cv, stop_cond = stop.cond * 10, mode_reg_A = mode_reg_A,...)
  # find a rough estimate of the attachment function based on the model Ak = k^alpha
  result_temp  <- PAFit(net_stat, 
                        mode_f      = "Log_linear"            , 
                        only_PA     = TRUE                    ,
                        alpha_start = cv_result$alpha_optimal , 
                        stop_cond   = stop.cond * 10          , # loose convergence condition 
                        ...) 
  # feed the estimated attachment function for a warm-start re-run with nonparametric attachment function
  result  <- PAFit(net_stat, 
                   r           = cv_result$r_optimal     , 
                   mode_f      = "Linear_PA"             ,  
                   only_PA     = TRUE                    ,
                   #alpha_start = result_temp$alpha   , 
                   alpha_start = cv_result$alpha_optimal , 
                   stop_cond   = stop.cond               ,
                   mode_reg_A  = mode_reg_A              ,  
                   ...) 
  
  return(list(cv_data = data_cv, cv_result = cv_result, estimate_result = result))
}