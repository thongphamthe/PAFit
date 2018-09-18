# function to summarize estimation results
summary.PAFit_result <- function(object,...){
    cat("Estimation results by the PAFit method. \n");
  
     if (object$only_PA == TRUE) {
        cat("Mode: Only the attachment function was estimated. \n")     
    } else if (object$only_f == TRUE) {
        cat("Mode: Only node fitnesses were estimated. \n")
    }
    else {
        cat("Mode: Both the attachment function and node fitness were estimated. \n")
    }
    #cat("Form of the PA function:",object$mode_f,"\n");
  
    if (object$only_f == FALSE) {
        if (object$auto_lambda == TRUE) {
            cat("Selected r parameter:", object$ratio,"\n");  
        } else cat("Lambda used:", object$lambda,"\n");
    }
    if (object$only_PA == FALSE)
        cat("Selected s parameter: ",object$shape,"\n", sep = "")
    cat("Estimated attachment exponent:",object$alpha,"\n");
    if (object$ci[1] == "N") {
      cat("No possible confidence interval for the estimated attachment exponent.\n");
    } else if (object$mode_f != "Log_linear") {
      cat("Attachment exponent ","\u00B1", " 2 s.d.", ": (", object$ci[1], ",", 
          object$ci[2],")\n",sep = "");
    }
    else {
      cat("Attachment exponent ","\u00B1", " 2 s.d.", ": (", object$ci[1], ",", 
          object$ci[2],")\n",sep = "");  
    }
    cat("-------------------------------------------\n")
    cat("Additional information: \n");
    cat("Number of bins:", object$g,"\n");
    cat("Number of iterations:",length(object$objective_value) - 1,"\n");
    cat("Stopping condition:", object$stop_cond,"\n");
}