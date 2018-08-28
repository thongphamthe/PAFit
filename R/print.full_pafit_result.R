# function to summarize estimation results
print.Full_PAFit_result <- function(x,...){
  object <- x
  cat("Estimation results by the PAFit method. \n");
  mode <- 0
  object_estimate_result <- object$estimate_result
  if (object_estimate_result$only_PA == TRUE) {
    cat("Mode: Only the attachment function was estimated. \n") 
    mode <- 0
  } else if (object_estimate_result$only_f == TRUE) {
    cat("Mode: Only node fitnesses were estimated. \n")
    mode <- 1
  }
  else {
    cat("Mode: Both the attachment function and node fitness were estimated. \n")
    mode <- 2
  }
  #cat("Form of the PA function:",object_estimate_result$mode_f,"\n");
  
  if (mode == 0 || mode == 2) {
    if (object_estimate_result$auto_lambda == TRUE) {
      cat("Selected r parameter:", object_estimate_result$ratio,"\n");
    } else cat("Lambda used:", object_estimate_result$lambda,"\n");
  }
  if (mode == 1 || mode == 2)
    cat("Selected s parameter:",object_estimate_result$shape,"\n")
  if (mode == 0 || mode == 2)
    cat("Estimated attachment exponent:",object_estimate_result$alpha,"\n");
  if (mode == 0 || mode == 2) {
    if (object_estimate_result$ci[1] == "N") {
      cat("No possible interval for the estimated attachment exponent.\n");
    } else if (object_estimate_result$mode_f != "Log_linear") {
      cat("Attachment exponent ","\u00B1", " 2 s.d.", ": (", object_estimate_result$ci[1], ",", 
          object_estimate_result$ci[2],")\n",sep = "");
    }
    else {
      cat("Attachment exponent ","\u00B1", " 2 s.d.", ": (", object_estimate_result$ci[1], ",", 
          object_estimate_result$ci[2],")\n",sep = "");  
    }
  }
}