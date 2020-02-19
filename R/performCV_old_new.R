.performCV_old_new <- function(cv_data                                            ,
                              r              = 10^c(-5 , -4 , -3 , -2, - 1, 0, 1) ,
                              s              = 10^c(- 1, 0 ,  1 ,  2 ,  3)        , 
                              stop_cond      = 10^-6                              ,
                              print_out      = FALSE                              ,
                              cv_deg_thresh  = c(1,10)                            ,
                              normal_start_f = TRUE                               ,
                              weight_f       = 0                                  ,
                              ...) { 
  
  num_deg_thresh <- length(cv_deg_thresh)
  # check possible
  for (dd in 1:length(cv_deg_thresh))
    while (TRUE) {
      chosen_node      <- names(cv_data$stats$z_j[cv_data$stats$z_j >= cv_deg_thresh[dd]])  
      if (length(chosen_node) <= 1)
        cv_deg_thresh[dd] <- round(cv_deg_thresh[dd]/2)   
      else break;
    }
  return(.performCV_core_new_new(cv_data,r = r, s = s, stop_cond = stop_cond, print_out = print_out,
                                 cv_deg_thresh = cv_deg_thresh, normal_start_f = normal_start_f,
                                 weight_f = weight_f))
}


