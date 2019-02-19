# function to generate simulated network  with PA and node fitness changes over time
.time_workflow <-
  function(num_period         = 50, # number of periods
           period_length      = 10,
           repeat_time        = 3 ,
           threshold          = 0.1,
           name_start = "Experiment_",
           directory  = "./") {
    options(scipen=999)
    #fit_list     <- vector(mode = "list", length = num_period)
    #var_fit_list <- vector(mode = "list", length = num_period)
    theta_list     <- vector(mode = "list", length = num_period)
    var_theta_list <- vector(mode = "list", length = num_period)
    k_list         <- vector(mode = "list", length = num_period)
    stats_list   <- vector(mode = "list", length = num_period)
    alpha_list   <- rep(0,num_period)
    upper_alpha  <- rep(0,num_period)
    lower_alpha  <- rep(0,num_period)
    r_list       <- rep(0,num_period)
    s_list       <- rep(0,num_period)
    z_list       <- vector(mode = "list",length = num_period) # the order of fitness
    load(file = paste0(directory,name_start,1,".Rdata"))
    n             <- length(node_id)
    appear        <- rep(0,n)
    names(appear) <- as.character(as.numeric(node_id))
    fit_matrix    <- matrix(-1,nrow = num_period, ncol = n)
    colnames(fit_matrix) <- as.character(as.numeric(node_id))
    var_f_matrix  <- fit_matrix
    
    for (ii in 1:num_period) {
        #print(ii)
        stats_part <- NULL
        net_part   <- NULL
        load(file = paste0(directory,name_start,ii,".Rdata"))
        
        stats_list[[ii]] <- stats_part
        result           <- joint_estimate(net_part,stats_part)
        as.character(as.numeric(stats_part$node_before_final))
        fit_matrix[ii,names(result$estimate_result$f)]   <- result$estimate_result$f
        var_f_matrix[ii,names(result$estimate_result$f)] <- result$estimate_result$var_f
        alpha_list[ii]  <- result$estimate_result$alpha
        upper_alpha[ii] <- result$estimate_result$ci[2]
        lower_alpha[ii] <- result$estimate_result$ci[1]
        r_list[ii]      <- result$estimate_result$ratio
        s_list[ii]      <- result$estimate_result$shape
         
        z_list[[ii]]       <- stats_part$z_j
        theta_list[[ii]]   <- result$estimate_result$theta
        var_theta_list[[ii]] <- result$estimate_result$var_theta
        k_list[[ii]]       <- result$estimate_result$center_k
    }
    node_id <- as.character(sort(as.numeric(unique(node_id))))
    if (repeat_time > 0) {
        for (t in 1:repeat_time){
        #print(paste0("Loop ",t))  
        # smoothing the fit matrix
            for (jj in 1:dim(fit_matrix)[2])
                fit_matrix[,jj] <- .smooth(series = fit_matrix[,jj],variance = var_f_matrix[,jj],
                                           threshold = threshold)  
            #print("After smoothing")  
        # then feed that smooth fitness to re-estimate PA
            for (ii in 1:num_period) {
            #load(file = paste0(directory,name_start,ii,".Rdata"))
            #save(stats_list,file = "temp.Rdata")
            result_A         <- PAFit(stats_list[[ii]], only_PA = "TRUE", 
                                      true_f = fit_matrix[ii,names(stats_list[[ii]]$z_j)],
                                       r = r_list[ii], s = s_list[ii])
            #save(result_A, file ="temp_A.Rdata")
            theta_list[[ii]]     <- result_A$theta
            alpha_list[ii]       <- result_A$alpha
            var_theta_list[[ii]] <- result_A$var_theta
            upper_alpha[ii]      <- result_A$ci[2]
            lower_alpha[ii]      <- result_A$ci[1]
            }
         # print("After re-estimate PA")  
        #then use the re-estimate PA to re-estimate fitness  
            for (ii in 1:num_period) {
            result_f           <- PAFit(stats_list[[ii]], only_f = "TRUE", true_A = theta_list[[ii]],
                                       r = r_list[ii], s = s_list[ii])
            fit_matrix[ii,names(result_f$f)]   <- result_f$f
            var_f_matrix[ii,names(result_f$f)] <- result_f$var_f
            }
        }
    }
    result <- list( fit_mat       = fit_matrix,
                    var_mat       = var_f_matrix,
                    theta_list    = theta_list,
                    center_k_list = k_list,
                    alpha_list    = alpha_list,
                    upper_alpha   = upper_alpha,
                    lower_alpha  = lower_alpha,
                    r_list      = r_list,
                    s_list      = s_list)
    return(result)
  }