.performCV_core_new_new <- function(cv_data                                       ,
                                   r              = 10^c(-5 , -4 , -3 , -2 , - 1) ,
                                   s              = 10^c(-1 ,  0 ,  1 ,  2 ,   3) , 
                                   stop_cond      = 10^-6                         ,
                                   print_out      = FALSE                               ,
                                   cv_deg_thresh  = c(1,10)                             ,
                                   normal_start_f = TRUE                                ,
                                   weight_f       = 0                                   ,
                                   ...) { 
  
  result_1 <- .one_cycle(cv_data, r, s,
                         estimated_fitness_start = NULL,
                         estimated_PA_start      = NULL, 
                         alpha_start             = NULL,
                         stop_cond,print_out,cv_deg_thresh,normal_start_f,weight_f);
  
  #plot(result_1$center_k,result_1$estimated_PA, log = "xy")
  #print(paste0("Max of theta outside one cycle: ",max(result_1$estimated_PA)))
  
  s <- c(0.5, 0.75, 1, 1.5, 2) * result_1$s_optimal;
  s <- s[s <= 10^4]
  s <- s[s >= 10^-2]
  r <- c(0.5, 0.75, 1, 1.5, 2) * result_1$r_optimal;
  r <- r[r >= 10^-5]
  r <- r[r <= 10^1]
  estimated_PA      <- result_1$estimated_PA;
  estimated_fitness <- result_1$estimated_fitness;
  estimated_alpha   <- result_1$alpha_optimal;
  
  #print(estimated_PA)
  #print(estimated_fitness)
  #print(estimated_alpha)
  
  result_2 <- .one_cycle(cv_data,r,s,estimated_fitness_start = estimated_fitness,
                         estimated_PA_start = estimated_PA, alpha_start = estimated_alpha, 
                         stop_cond, print_out, cv_deg_thresh,
                         normal_start_f, weight_f);
  
  
  
  #exit();
  
  s <- c(1/1.75 , 1/1.5 , 1/1.25 ,1, 1.25 , 1.5 , 1.75)  * result_2$s_optimal;
  s <- s[s <= 10^4]
  s <- s[s >= 10^-2]
  
  r <-  c(1/1.75 , 1/1.5 , 1/1.25 ,1, 1.25 , 1.5 , 1.75) * result_2$r_optimal;
  r <- r[r >= 10^-5]
  r <- r[r <= 10^1]
  
  estimated_PA      <- result_2$estimated_PA;
  estimated_fitness <- result_2$estimated_fitness;
  estimated_alpha   <- result_2$alpha_optimal;
  
 
  
  result_2_5 <- .one_cycle(cv_data,r,s,
                           estimated_fitness_start = estimated_fitness,
                           estimated_PA_start = estimated_PA, alpha_start = estimated_alpha, 
                           stop_cond, print_out, cv_deg_thresh,
                           normal_start_f, weight_f);
  
  s <- c(1/1.75 , 1/1.5 , 1/1.25 ,1, 1.25 , 1.5 , 1.75)  * result_2_5$s_optimal;
  s <- s[s <= 10^4]
  s <- s[s >= 10^-2]
  
  r <-  c(1/1.75 , 1/1.5 , 1/1.25 ,1, 1.25 , 1.5 , 1.75) * result_2_5$r_optimal;
  r <- r[r >= 10^-5]
  r <- r[r <= 10^1]
  
  result_2_75 <- .one_cycle(cv_data,r,s,
                           estimated_fitness_start = estimated_fitness,
                           estimated_PA_start = estimated_PA, alpha_start = estimated_alpha, 
                           stop_cond, print_out, cv_deg_thresh,
                           normal_start_f, weight_f);
  
 
  
  
  s <- c(1/1.75 , 1/1.5 , 1/1.25 ,1, 1.25 , 1.5 , 1.75)  * result_2_75$s_optimal;
  s <- s[s <= 10^4]
  s <- s[s >= 10^-2]
  
  r <-  c(1/1.75 , 1/1.5 , 1/1.25 ,1, 1.25 , 1.5 , 1.75) * result_2_75$r_optimal;
  r <- r[r >= 10^-5]
  r <- r[r <= 10^1]
  
  estimated_PA      <- result_2_75$estimated_PA;
  estimated_fitness <- result_2_75$estimated_fitness;
  estimated_alpha   <- result_2_75$alpha_optimal;
  
  result_3 <- .one_cycle(cv_data,r,s,
                         estimated_fitness_start = estimated_fitness,
                         estimated_PA_start = estimated_PA, alpha_start = estimated_alpha, 
                         stop_cond, print_out, cv_deg_thresh,
                         normal_start_f, weight_f);
  
  r_optimal         <- result_3$r_optimal
  lambda_optimal    <- result_3$lambda_optimal
  s_optimal         <- result_3$s_optimal
  alpha_optimal     <- result_3$alpha_optimal
  estimated_fitness <- result_3$estimated_fitness
  estimated_PA      <- result_3$estimated_PA
  
  
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


