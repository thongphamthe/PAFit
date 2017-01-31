performCV <- function(cv_data,
                      r              = 10^c(-2,-1,0,1,2),
                      s              = 10^c(-1,1,2,3,4) , 
                      stop_cond      = 10^-7            , 
                      only_PAFit     = TRUE             , 
                      silent         = FALSE            , 
                      only_loglinear = FALSE            ,
                      ...) { 
  ratio_vec_PAFit        <- r
  rate_PAFit             <- s
  ratio_vec_PA           <- ratio_vec_PAFit
  rate_Fit               <- rate_PAFit  
  PAFit_each             <- matrix(0,nrow = length(ratio_vec_PAFit), ncol = length(rate_PAFit))
  rownames(PAFit_each)   <- ratio_vec_PAFit
  colnames(PAFit_each)   <- rate_PAFit
  Fit_each_linear        <- rep(0,length(rate_Fit))
  names(Fit_each_linear) <- rate_Fit
  PA_each                <- rep(0,length(ratio_vec_PA))
  names(PA_each)         <- ratio_vec_PA
  Fit_each               <- rep(0,length(rate_Fit))
  names(Fit_each)        <- rate_Fit
  alpha_each             <- rep(0,length(rate_PAFit))
  names(alpha_each)      <- rate_PAFit
  FitMultinomial         <- function(true,dat){
      true[true == 0] <- 1
      return(sum(dat*log(true)))
  }
  count <- 0
  total <- length(ratio_vec_PAFit) * length(rate_PAFit) + (!only_PAFit) *
           (2 * length(rate_PAFit) + length(ratio_vec_PAFit) + (!only_loglinear) * length(rate_PAFit))[1]
  #Full model
  if (!only_loglinear)
  for (i in 1:length(ratio_vec_PAFit))
    for (j in 1:length(rate_PAFit)){
      count <- count + 1
      if (silent == FALSE)
          print(paste0("Processing case ",count, " of ",total))
      result_PAFit <- PAFit(cv_data$stats,s = rate_PAFit[j], r = ratio_vec_PAFit[i], auto_stop =  TRUE, 
                            stop_cond = stop_cond, normalized_f = FALSE,...)   
      for (k in 1:length(cv_data$m_each))
        if (cv_data$m_each[k] != 0) { 
          prob_PAFit      <- result_PAFit$A[cv_data$deg_each[k,] + 1]* result_PAFit$f[as.character(cv_data$stats$f_position)]; 
          prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE); 
          prob_PAFit[sapply(prob_PAFit,is.na)] <- 0; 
          PAFit_each[i,j] <-PAFit_each[i,j] + 
            FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,] * cv_data$m_each[k])) ;
        }
    }
  #linear A
  if (!only_loglinear)
  if (FALSE == only_PAFit)
  if (length(rate_Fit) > 0) {
    for (i in 1:length(rate_Fit)) {
      count <- count + 1
      if (silent == FALSE)
          print(paste0("Processing case ",count, " of ",total))
      result_PAFit <- PAFit(cv_data$stats, mode_f = "Linear_PA", only_f = TRUE,s = rate_Fit[i], auto_stop =  TRUE, 
                            stop_cond = stop_cond,normalized_f = FALSE,...)     
      for (k in 1:length(cv_data$m_each)) 
        if (cv_data$m_each[k] != 0) {
          prob_PAFit      <- result_PAFit$f[as.character(cv_data$stats$f_position)] *
                             pmax(cv_data$deg_each[k,],1)  
          prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE)
          prob_PAFit[sapply(prob_PAFit,is.na)] <- 0
          Fit_each_linear[i] <-Fit_each_linear[i] + FitMultinomial(true = as.vector(prob_PAFit), 
                                                                   dat = as.vector(cv_data$prob_em_each[k,] * cv_data$m_each[k]))
        }
      
    }  
  }
  if (!only_loglinear)
  #Only PA
  if (FALSE == only_PAFit)
  if (length(ratio_vec_PA) > 0) {
    for (i in 1:length(ratio_vec_PA)) {
      count <- count + 1
      if (silent == FALSE)
          print(paste0("Processing case ",count, " of ",total))
      result_PAFit <- PAFit(cv_data$stats, only_PA = TRUE, auto_lambda = TRUE, r = ratio_vec_PA[i], auto_stop =  TRUE, 
                            stop_cond = stop_cond,normalized_f = FALSE,...)
      
      for (k in 1:length(cv_data$m_each)) 
        if (cv_data$m_each[k] != 0) {
          prob_PAFit      <- result_PAFit$A[cv_data$deg_each[k,] + 1]
          prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE)
          prob_PAFit[sapply(prob_PAFit,is.na)] <- 0
          PA_each[i] <- PA_each[i] + FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,] * cv_data$m_each[k]))
        }
    }
  }
  # Only f
  if (!only_loglinear)
  if (FALSE == only_PAFit)
  if (length(rate_Fit) > 0) {
    for (i in 1:length(rate_Fit)) {
      count <- count + 1
      if (silent == FALSE)
          print(paste0("Processing case ",count, " of ",total))
      result_PAFit <- PAFit(cv_data$stats, mode_f = "Constant_PA", only_f = TRUE,s = rate_Fit[i], 
                            auto_stop =  TRUE, 
                            stop_cond = stop_cond,normalized_f = FALSE,...)     
      for (k in 1:length(cv_data$m_each)) 
        if (cv_data$m_each[k] != 0) {
          prob_PAFit      <- result_PAFit$f[as.character(cv_data$stats$f_position)] 
          prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE)
          prob_PAFit[sapply(prob_PAFit,is.na)] <- 0
          Fit_each[i]     <-Fit_each[i] + FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,] * cv_data$m_each[k]))
        }
    }
  }
  #k^alpha
  if (only_loglinear)
  if (FALSE == only_PAFit)
  for (j in 1:length(rate_PAFit)) {
    count <- count + 1
    if (silent == FALSE)
        print(paste0("Processing case ",count, " of ",total))
    result_PAFit <- PAFit(cv_data$stats, mode_f = "Log_linear",s = rate_PAFit[j],
                          auto_stop =  TRUE, stop_cond = stop_cond,normalized_f = FALSE,...);
    for (k in 1:length(cv_data$m_each)) 
      if (cv_data$m_each[k] != 0) {
       
        
        prob_PAFit      <- result_PAFit$A[cv_data$deg_each[k,] + 1] * result_PAFit$f[as.character(cv_data$stats$f_position)] 
        
        prob_PAFit      <- prob_PAFit / sum(prob_PAFit,na.rm = TRUE) 
        
        prob_PAFit[sapply(prob_PAFit,is.na)] <- 0 
        
        alpha_each[j] <- alpha_each[j] + 
          FitMultinomial(true = as.vector(prob_PAFit), dat = as.vector(cv_data$prob_em_each[k,] * cv_data$m_each[k])) 
      }
  }
  max_id    <- which.max(PAFit_each)[1]
  r_index   <- max_id %% length(r)
  if (r_index == 0)
      r_index <-  1
  s_index   <- ceiling(max_id / length(s))
  r_optimal <- r[r_index]
  #print(r)
  #print(r_index)
  #print(r_optimal)
  s_optimal <- s[s_index]
  result    <- list(PAFit_each      = PAFit_each, 
                    Fit_each_linear = Fit_each_linear, 
                    PA_each         = PA_each, 
                    Fit_each        = Fit_each, 
                    alpha_each      = alpha_each, 
                    r_optimal       = r_optimal, 
                    s_optimal       = s_optimal)
  class(result) <- "CV_Result"
  if(FALSE == silent) {
      print(paste0("Optimal r parameter is: ",r_optimal));
      print(paste0("Optimal s parameter is: ",s_optimal));
  }
  return(result)
}


