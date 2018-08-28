summary.Linear_PA_test_result <- function(object,...){
  cat("\nContaining the AIC and BIC of fitting various distributions to the degree distribution.");
  full_result             <-  object
  chosen_model_AIC        <-  vector(mode = "list", length = 5)
  chosen_model_BIC        <-  vector(mode = "list", length = 5)
  names(chosen_model_AIC) <- c("yule","waring","nb","geom","pois")
  names(chosen_model_BIC) <- c("yule","waring","nb","geom","pois")
  AIC_vector              <- 1:5
  BIC_vector              <- 1:5
  k_min                    <- 1:length(full_result[[1]])
  for (j in 1:5) {
    min_AIC_index <- 1
    min_BIC_index <- 1
    if (length(k_min) > 1)
      for (i in 2:length(k_min)) {
        if (full_result[[j]][[i]]$AIC < full_result[[j]][[min_AIC_index]]$AIC)
          min_AIC_index <- i  
        if (full_result[[j]][[i]]$BIC < full_result[[j]][[min_BIC_index]]$BIC)
          min_BIC_index <- i  
      }
    chosen_model_AIC[[j]] <- full_result[[j]][[min_AIC_index]]
    chosen_model_BIC[[j]] <- full_result[[j]][[min_BIC_index]]
    AIC_vector[j] <- full_result[[j]][[min_AIC_index]]$AIC
    BIC_vector[j] <- full_result[[j]][[min_BIC_index]]$BIC
  }
  data_frame <- matrix(0,nrow = 5, ncol = 6)
  for (i in 1:length(order(AIC_vector))) {
    data_frame[i,1] <- names(chosen_model_AIC)[order(AIC_vector)[i]]
    data_frame[i,2] <- chosen_model_AIC[[order(AIC_vector)[i]]]$k_min
    data_frame[i,3] <- chosen_model_AIC[[order(AIC_vector)[i]]]$dim
    data_frame[i,4] <- format(round(chosen_model_AIC[[order(AIC_vector)[i]]]$val,2),nsmal = 2)
    data_frame[i,5] <- format(round(chosen_model_AIC[[order(AIC_vector)[i]]]$AIC,2),nsmal = 2)
    data_frame[i,6] <- format(round(chosen_model_AIC[[order(AIC_vector)[i]]]$BIC,2),nsmal = 2)
  }
  colnames(data_frame) <-  c("Model","k_min","Dimension","Log-Likelihood","AIC","BIC")
  
  cat("\nTop models ordered by the value of the AIC:");
  temp <- knitr::kable(head(data_frame))
  print(temp)
  # cat("\nTop five models ordered by the value of the BIC:");
  # data_frame <- matrix(0,nrow = 5, ncol = 6)
  # for (i in 1:length(order(BIC_vector))) {
  #   data_frame[i,1] <- names(chosen_model_BIC)[order(BIC_vector)[i]]
  #   data_frame[i,2] <- chosen_model_BIC[[order(BIC_vector)[i]]]$k_min
  #   data_frame[i,3] <- chosen_model_BIC[[order(BIC_vector)[i]]]$dim
  #   data_frame[i,4] <- chosen_model_BIC[[order(BIC_vector)[i]]]$val
  #   data_frame[i,5] <- chosen_model_BIC[[order(BIC_vector)[i]]]$AIC
  #   data_frame[i,6] <- chosen_model_BIC[[order(BIC_vector)[i]]]$BIC
  # }
  # colnames(data_frame) <-  c("Model","k_min","Dimension","Log-Likelihood","AIC","BIC")
  # temp <- knitr::kable(head(data_frame))
  # print(temp)
}