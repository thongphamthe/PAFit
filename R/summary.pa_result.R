summary.PA_result <- function(object,...){
  cat("\nPA_result object contains estimated attachment function. \n");
  cat("Number of bins: ", object$G,"\n");
  cat("Estimated attachment exponent:", object$alpha, "\n");
  
  if (!is.null(object$stop_cond)) {
      cat("Stopping condition:", object$stop_cond,"\n");
      cat("Auto Lambda: ",object$auto_lambda,"\n");
      if (object$auto_lambda == TRUE) 
          cat("Ratio (r): ", object$ratio,"\n");
      cat("Lambda used: ", object$lambda,"\n");
  }
}