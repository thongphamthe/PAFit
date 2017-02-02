# function to summarize estimation results
summary.PAFit_result <- function(object,...){
    cat("\nPAFit_result object \n");
    cat("Estimation results by the PAFit method. \n");
    cat("Mode:","Only PA: ",object$only_PA,"; Only f:", object$only_f,"\n");
    cat("Form of the PA function:",object$mode_f,"\n");
    if (object$auto_lambda == TRUE) {
        cat("Ratio (r): ", object$ratio,"\n");
    } else cat("Lambda used: ", object$lambda,"\n");
    cat("Prior of node fitness: shape: ",object$shape,"; rate: ",object$rate,"\n")
    cat("Estimated attachment exponent: ",object$alpha,"\n");
    if (object$ci[1] == "N") {
      cat("No possible confidence interval for the estimated attachment exponent.\n");
    } else if (object$mode_f != "Log_linear") {
      cat("95% confidence interval of the attachment exponent: (", object$ci[1], ",", 
          object$ci[2],")\n");
    }
    else {
      cat("Two-sigma confidence interval of the attachment exponent: (", object$ci[1], ",", 
          object$ci[2],")\n");  
    }
    cat("-------------------------------------------\n")
    cat("Additional information: \n");
    cat("Number of bins: ", object$G,"\n");
    cat("Threshold of the number of new edges acquired: ", object$deg_thresh,"\n");
    cat("Number of nodes satisfied the threshold: ",sum(object$f != 1),"\n");
    cat("Number of iterations: ",length(object$objective_value) - 1,"\n");
    cat("Stopping condition:", object$stop_cond,"\n");
}