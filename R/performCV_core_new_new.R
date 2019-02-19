.performCV_core_new_new <- function(cv_data                                     ,
                                   r                   ,
                                   s            , 
                              
                                   stop_cond      = 10^-6                       ,
                                   print_out      = FALSE                               ,
                                   cv_deg_thresh  = c(1,10)                             ,
                                   normal_start_f = TRUE                                ,
                                   weight_f       = 0                                   ,
                                   ...) { 

  multiply_vector_s <- c(1/4,1/3,1/2,1,2,3,4) 
  multiply_vector_r <- c(1/4,1/3,1/2,1,2,3,4)
  r_optima_old <- NULL
  s_optima_old <- NULL
  break_flag   <- FALSE
  count        <- 0
  for (i in 1:20) {
       if (i == 1) {
       result <- .one_cycle(cv_data, r, s,
                         estimated_fitness_start = NULL,
                         estimated_PA_start      = NULL, 
                         alpha_start             = NULL,
                         stop_cond,print_out,cv_deg_thresh,normal_start_f,weight_f);
       } else {
      result <- .one_cycle(cv_data,r,s,estimated_fitness_start = estimated_fitness,
                    estimated_PA_start = estimated_PA, alpha_start = estimated_alpha, 
                    stop_cond, print_out, cv_deg_thresh,
                    normal_start_f, weight_f);   
       }
      if (i > 1) {
          if ((s_optima_old == result$s_optimal) & (r_optima_old == result$r_optimal))  {
              if (count < 2) {  
                  multiply_vector_s <- c(0.75,0.85,0.95,1,1.25,1.5,1.75)
                  multiply_vector_r <- c(0.75,0.85,0.95,1,1.25,1.5,1.75)
                  count             <- count + 1
              } else {
                break_flag <- TRUE  
              }
          }
      }
    estimated_PA      <- result$estimated_PA;
    estimated_fitness <- result$estimated_fitness;
    estimated_alpha   <- result$alpha_optimal;
     if (break_flag) break;
      r_optima_old <- result$r_optimal
      s_optima_old <- result$s_optimal
      s <- multiply_vector_s * result$s_optimal;
      s <- s[s <= 10^3]
      s <- s[s > 1]
      r <- multiply_vector_r * result$r_optimal;
      r <- r[r >= 10^-5]
      r <- r[r <= 10]
  }
  
  r_optimal         <- result$r_optimal
  lambda_optimal    <- result$lambda_optimal
  s_optimal         <- result$s_optimal
  alpha_optimal     <- result$alpha_optimal
  estimated_fitness <- result$estimated_fitness
  estimated_PA      <- result$estimated_PA
  
  
  result    <- list(r_optimal         = r_optimal, 
                    lambda_optimal    = lambda_optimal, 
                    s_optimal         = s_optimal,
                    alpha_optimal     = alpha_optimal,
                    estimated_fitness = estimated_fitness,
                    estimated_PA      = estimated_PA,
                    cv_deg_thresh     = cv_deg_thresh)
  
  class(result) <- "CV_Result"
  return(result)
}


