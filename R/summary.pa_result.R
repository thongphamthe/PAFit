summary.PA_result <- function(object,...){
  cat("\nPA_result object contains estimated attachment function. \n");
  cat("Number of bins: ", object$G,"\n");
  cat("Estimated attachment exponent:", object$alpha, "\n");
  cat("95% confidence interval of the attachment exponent: (", object$ci[1], ",", 
      object$ci[2],")\n");
}