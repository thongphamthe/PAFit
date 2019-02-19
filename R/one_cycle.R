.one_cycle <- function(cv_data                            ,
                       r                                  ,
                       s                                  ,
                       estimated_fitness_start = NULL     ,
                       estimated_PA_start      = NULL     ,
                       alpha_start             = NULL     ,
                       stop_cond               = 10^-6    ,
                       print_out               = FALSE    ,
                       cv_deg_thresh           = c(1,10)  ,
                       normal_start_f          = TRUE     ,
                       weight_f                = 0        ,
                       ...) { 
  
  FitMultinomial         <- function(true,dat){
    temp <- dat*log10(true); 
    if (sum(is.na(temp)) > 0) {
      #temp[is.na(temp)] <- 0;
      return(-Inf);
    } else { return(sum(temp));}
  }
  rate_PAFit             <- s
  
  if (is.null(estimated_PA_start)) {
      alpha_optimal        <- rep(-1,length(cv_deg_thresh))
  }
  
  alpha_each           <- matrix(0,nrow = length(cv_deg_thresh), ncol = length(rate_PAFit))
  colnames(alpha_each) <- rate_PAFit
  estimated_fitness    <- list(length = length(cv_deg_thresh))
  #estimated_PA        <- NULL
  max_val              <- rep(-Inf,length(cv_deg_thresh))
  s_optimal            <- rep(0,length(cv_deg_thresh))
  
  
  #print(rate_PAFit)
  
  for (j in 1:length(rate_PAFit)) {
    
    #print(paste0("Processing case ",count, " of ",total))
    if (is.null(estimated_PA_start)) {
        if (!is.null(estimated_fitness_start)) { 
            result_PAFit  <- PAFit(cv_data$stats , mode_f    = "Log_linear"    , s = rate_PAFit[j], 
                                   auto_stop      =  TRUE , stop_cond  = stop_cond , normalized_f = FALSE ,
                                   normal_start_f = normal_start_f,
                                   weight_power   = weight_f, start_f = estimated_fitness_start,
                                   ...);
        }
        else 
            result_PAFit  <- PAFit(cv_data$stats     , mode_f    = "Log_linear"    , s = rate_PAFit[j], 
                                   auto_stop      =  TRUE , stop_cond  = stop_cond , normalized_f = FALSE ,
                                   normal_start_f = normal_start_f,
                                   weight_power   = weight_f, 
                                    ...);  
    }
    else {
        #print("inside here");  
        if (!is.null(estimated_fitness)) { 
            result_PAFit  <- PAFit(cv_data$stats     , true_A = estimated_PA_start, 
                                   s = rate_PAFit[j] , only_f = TRUE,
                                   auto_stop =  TRUE , 
                                   stop_cond  = stop_cond , normalized_f = FALSE ,
                                   normal_start_f = normal_start_f,
                                   weight_power   = weight_f, start_f = estimated_fitness_start,
                                   ...);  
            #plot(result_PAFit$center_k,result_PAFit$theta,log = "xy", main = "Inside loop");
            #print(paste0("Max of theta inside one cycle: ",max(estimated_PA_start)))
        }
        else 
            result_PAFit  <- PAFit(cv_data$stats     , true_A = estimated_PA_start, 
                                   s = rate_PAFit[j] , only_f = TRUE,
                                   auto_stop =  TRUE , 
                                   stop_cond  = stop_cond , normalized_f = FALSE ,
                                   normal_start_f = normal_start_f,
                                   weight_power   = weight_f,
                                   ...);    
    }
    if (is.null(estimated_PA_start)) {
        alpha_temp    <- result_PAFit$alpha
        #print(alpha_temp)
        if (TRUE == result_PAFit$diverge_zero) {
            #print("Diverge to zero");
            alpha_each[,j] <- -Inf
        }
        else  { 
            for (dd in 1:length(cv_deg_thresh)) { 
                alpha_each[dd,j]  <- 0 
                name_before_final <- intersect(as.character(as.numeric(cv_data$stats$node_before_final)),names(cv_data$stats$z_j))
                z_j_before_final  <- cv_data$stats$z_j[name_before_final]
                chosen_node_big   <- names(z_j_before_final[z_j_before_final >= cv_deg_thresh[dd]])
                #print(colnames(cv_data$deg_each))
                #print(chosen_node_big)
            for (k in 1:length(cv_data$m_each)) 
                if (cv_data$m_each[k] != 0) {
                #PA              <- result_PAFit$A[cv_data$deg_each[k,] + 1] 
                #PA[is.na(PA)]   <- result_PAFit$A[length(result_PAFit$A)]
                chosen_node     <- chosen_node_big[cv_data$deg_each[k,chosen_node_big] != 0]   
                PA              <- pmax(cv_data$deg_each[k,chosen_node],1)^alpha_temp
                #print(PA)
                if (sum(PA == 0) > 0)
                    stop("sum(PA == 0) > 0")  
                #PA[PA == 0]     <- 1
                fitness         <- rep(0,dim(cv_data$deg_each[,chosen_node, drop = FALSE])[2])
                names(fitness)  <- colnames(cv_data$deg_each[,chosen_node, drop = FALSE])
                fitness[chosen_node] <- result_PAFit$f[chosen_node] 
               
                non_na_chosen   <- chosen_node[!is.na(PA) &  !is.na(fitness)]
                if (sum(is.na(PA)) > 0)
                  print(paste0("NA in PA: ",sum(is.na(PA)))) 
                if (sum(is.na(fitness)) > 0)
                  print(paste0("NA in fitness: ",sum(is.na(fitness))))
                prob_PAFit      <- PA[non_na_chosen] * fitness[non_na_chosen]
                prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
                #prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
                alpha_each[dd,j]   <- alpha_each[dd, j] + 
                                      FitMultinomial(true = as.vector(prob_PAFit), 
                                                     dat = as.vector(cv_data$prob_em_each[k,non_na_chosen])) 
            }
        }
      }
        for (dd in 1:length(cv_deg_thresh)) {
          #print(alpha_each[dd,j]);
          if (alpha_each[dd,j] > max_val[dd]) {
            estimated_fitness[[dd]]  <- result_PAFit$f[as.character(as.numeric(cv_data$stats$f_position))]
            s_optimal[dd]            <- rate_PAFit[j]
            max_val[dd]              <- alpha_each[dd,j]
            alpha_optimal[dd]        <- result_PAFit$alpha
          }
        }
    } else {
         ### estimated_PA is provided
          if (TRUE == result_PAFit$diverge_zero) {
              alpha_each[,j] <- -Inf
          } else  { 
             
              for (dd in 1:length(cv_deg_thresh)) { 
                  alpha_each[dd,j] <- 0  
                  name_before_final <- intersect(as.character(as.numeric(cv_data$stats$node_before_final)),names(cv_data$stats$z_j))
                  z_j_before_final  <- cv_data$stats$z_j[name_before_final]
                  chosen_node_big  <- names(z_j_before_final[z_j_before_final >= cv_deg_thresh[dd]])
                  for (k in 1:length(cv_data$m_each)) 
                  if (cv_data$m_each[k] != 0) {
                      chosen_node_ok_flag <- 1
                      chosen_node      <- chosen_node_big[!is.na(result_PAFit$A[as.character(cv_data$deg_each[k,chosen_node_big] + 1)])]
                      if (length(chosen_node) > 0) {
                          chosen_node      <- chosen_node[result_PAFit$A[as.character(cv_data$deg_each[k,chosen_node] + 1)] != 0]
                          if (length(chosen_node) > 0) {
                
                              PA              <- pmax(cv_data$deg_each[k,chosen_node],1)^alpha_start
                              if (sum(is.na(PA)) > 0) {
                              #PA[is.na(PA)] <- max(estimated_PA_start)  
                              #PA[is.na(PA)]   <- cv_data$deg_each[k,chosen_node][is.na(PA)]^alpha_start
                          
                                 stop("sum(is.na(PA)) > 0")
                             }
                             if (sum(PA == 0) > 0)
                                PA[PA == 0]     <- 1
                          } else {
                              chosen_node_ok_flag <- 0  
                          }
                      } else {
                        chosen_node_ok_flag <- 0  
                      }
                      if (chosen_node_ok_flag == 0) {
                          chosen_node <- chosen_node_big
                          PA          <- pmax(cv_data$deg_each[k,chosen_node],1)^alpha_start
                      }
                      fitness         <- rep(0,dim(cv_data$deg_each[,chosen_node, drop = FALSE])[2])
                      names(fitness)  <- colnames(cv_data$deg_each[,chosen_node, drop = FALSE])
                      fitness[chosen_node] <- result_PAFit$f[chosen_node] 
              
                      non_na_chosen   <- chosen_node[!is.na(PA) &  !is.na(fitness)]
                      if (sum(is.na(PA)) > 0)
                        print(paste0("NA in PA: ",sum(is.na(PA)))) 
                      if (sum(is.na(fitness)) > 0)
                        print(paste0("NA in fitness: ",sum(is.na(fitness))))
                      
                      prob_PAFit      <- PA[non_na_chosen] * fitness[non_na_chosen]
                      prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
                      #prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
                      alpha_each[dd,j]   <- alpha_each[dd, j] + 
                        FitMultinomial(true = as.vector(prob_PAFit), 
                                       dat = as.vector(cv_data$prob_em_each[k,non_na_chosen])) 
                 }
            }
          }
          for (dd in 1:length(cv_deg_thresh)) 
              if (alpha_each[dd,j] > max_val[dd]) {
                  estimated_fitness[[dd]]  <- result_PAFit$f[as.character(as.numeric(cv_data$stats$f_position))]
                  s_optimal[dd]            <- rate_PAFit[j]
                  max_val[dd]              <- alpha_each[dd,j]
                 
              }
          }
  }
  
  #print(alpha_each)
  
  #print(alpha_optimal)
  #print(alpha_each)
  
  ok_index                <- 1:length(cv_deg_thresh)
  
  s_optimal_final         <- mean(s_optimal[ok_index])
 
  #if (s_optimal_final == 0)
  
  #s_optimal_final         <- 10^(1/length(s_optimal[ok_index])* log10(prod(s_optimal[ok_index])))
  
  if (is.null(estimated_PA_start))
      alpha_optimal_final     <- mean(alpha_optimal[ok_index])
  estimated_fitness_final <- estimated_fitness[[ok_index[1]]]
  
  if (length(ok_index) > 1)
    for (temp_jj in 2:length(ok_index))
      estimated_fitness_final <- estimated_fitness_final + estimated_fitness[[ok_index[temp_jj]]] 
  estimated_fitness_final <- estimated_fitness_final / length(ok_index)
  
  #print("Reach here")

  ### One pass to find optimal r #######

  name_before_final <- intersect(as.character(as.numeric(cv_data$stats$node_before_final)),names(cv_data$stats$z_j))
  z_j_before_final  <- cv_data$stats$z_j[name_before_final]
  chosen_node_big  <- names(z_j_before_final[z_j_before_final >= cv_deg_thresh[1]])

  ratio_vec_PAFit        <- sort(r,decreasing = TRUE)
  PA_each                <- rep(0,length(ratio_vec_PAFit))
  names(PA_each)         <- ratio_vec_PAFit
  max_val                <- -Inf 
  
  
  estimated_PA           <- NULL
  
  switch_flag            <- 0
  
  for (i in 1:length(ratio_vec_PAFit)) {
  # print(i)
  #print(!is.null(estimated_PA))
  #print(is.null(estimated_PA_start))
    
    #print(paste0("Processing case ",count, " of a maximum of ",total))
    if (!is.null(estimated_PA))
      result_PAFit <- PAFit(cv_data$stats, 
                            s = s_optimal_final, 
                            r = ratio_vec_PAFit[i], 
                            only_PA   = TRUE, 
                            true_f    = estimated_fitness_final, 
                            auto_stop = TRUE, 
                            start_A   = estimated_PA,
                            stop_cond = stop_cond, 
                            normal_start_f = normal_start_f,
                            weight_power   = weight_f ,
                            normalized_f = FALSE, ...)
    
    else {
        if (is.null(estimated_PA_start)) {
            ## error here
            result_PAFit <- PAFit(cv_data$stats, 
                                  s = s_optimal_final, 
                                  r = ratio_vec_PAFit[i], 
                                  only_PA = TRUE, 
                                  true_f = estimated_fitness_final, 
                                  auto_stop =  TRUE, 
                                  alpha_start = alpha_optimal_final,
                                  stop_cond = stop_cond, 
                                  normal_start_f = normal_start_f,
                                  normalized_f = FALSE, ...)
        }
        else 
          result_PAFit <- PAFit(cv_data$stats, 
                                s = s_optimal_final, 
                                r = ratio_vec_PAFit[i], 
                                only_PA = TRUE, 
                                true_f = estimated_fitness_final, 
                                auto_stop =  TRUE, 
                                start_A = estimated_PA_start,
                                stop_cond = stop_cond, 
                                normal_start_f = normal_start_f,
                                normalized_f = FALSE, ...)  
    }
    
    
    alpha_temp    <- result_PAFit$alpha  
    
    
    if (TRUE == result_PAFit$diverge_zero)
        PA_each[i] <- -Inf
    else for (k in 1:length(cv_data$m_each))
        if (cv_data$m_each[k] != 0) { 
            chosen_node_ok_flag <- 1
            chosen_node    <- chosen_node_big[!is.na(result_PAFit$A[as.character(cv_data$deg_each[k,chosen_node_big] + 1)])]
            if (length(chosen_node) > 0) {
                chosen_node    <- chosen_node[result_PAFit$A[as.character(cv_data$deg_each[k,chosen_node] + 1)] != 0]
                if (length(chosen_node) > 0) {
                    PA             <- pmax(cv_data$deg_each[k,chosen_node],1)^alpha_temp
                    #PA             <- result_PAFit$A[as.character(cv_data$deg_each[k,chosen_node] + 1)]
                    if (sum(is.na(PA)) > 0) {
                    #print("it should not happen here")  
                        PA[is.na(PA)]   <- pmax(cv_data$deg_each[k,chosen_node][is.na(PA)],1)^alpha_temp
                    }
                    if (sum(PA == 0) > 0) {
                        #print("This should not happen")  
                       PA[PA == 0]     <- pmax(cv_data$deg_each[k,chosen_node][PA == 0],1)^alpha_temp
                    }
                } else {
                  chosen_node_ok_flag <- 0  
                }
            } else {
              chosen_node_ok_flag <- 0      
            }
            if (chosen_node_ok_flag == 0) {
               chosen_node <- chosen_node_big  
               PA          <- pmax(cv_data$deg_each[k,chosen_node],1)^alpha_temp
            }
            
            fitness         <- rep(0,dim(cv_data$deg_each[,chosen_node, drop = FALSE])[2])
            names(fitness)  <- colnames(cv_data$deg_each[,chosen_node, drop = FALSE])
            fitness[chosen_node] <- result_PAFit$f[chosen_node] 
        
            non_na_chosen   <- chosen_node[!is.na(PA) &  !is.na(fitness)]
            
            if (sum(is.na(PA)) > 0)
                print(paste0("NA in PA: ",sum(is.na(PA)))) 
            if (sum(is.na(fitness)) > 0)
                print(paste0("NA in fitness: ",sum(is.na(fitness))))
            
            prob_PAFit      <- PA[non_na_chosen] * fitness[non_na_chosen]
            prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
            #prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
            PA_each[i]      <- PA_each[i] + 
              FitMultinomial(true = as.vector(prob_PAFit), 
                             dat = as.vector(cv_data$prob_em_each[k,non_na_chosen])) 
        
         }
    if (i == 1) {  
         lambda_optimal    <- result_PAFit$lambda  
         r_optimal         <- ratio_vec_PAFit[i]
         max_val           <- PA_each[i]
         estimated_PA      <- result_PAFit$theta
         alpha_op_final    <- result_PAFit$alpha
    } else if (PA_each[i] > max_val) {
          lambda_optimal    <- result_PAFit$lambda  
          r_optimal         <- ratio_vec_PAFit[i]
          max_val           <- PA_each[i]
          estimated_PA      <- result_PAFit$theta  
          alpha_op_final    <- result_PAFit$alpha
    }
  }
  #if (print_out == TRUE)
  
  #print(paste0("s_optimal = ",s_optimal_final, " ,r_optimal = ", r_optimal)) 
  #print(paste0(s_optimal_final,";",r_optimal))
  
  result    <- list(r_optimal         = r_optimal, 
                    lambda_optimal    = lambda_optimal, 
                    s_optimal         = s_optimal_final,
                    estimated_fitness = estimated_fitness_final,
                    estimated_PA      = estimated_PA,
                    cv_deg_thresh     = cv_deg_thresh,
                    alpha_optimal     = alpha_op_final,
                    
                    center_k          = result_PAFit$center_k)
  class(result) <- "CV_Result"
  return(result)
}