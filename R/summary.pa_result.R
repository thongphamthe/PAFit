summary.PA_result <- function(object,...){
  cat("\nContaining the estimated attachment function. \n");
  cat("Number of bins: ", object$g,"\n");
  cat("Estimated attachment exponent:", object$alpha, "\n");
  if (object$ci[1] == "N") {
      cat("No possible confidence interval for the estimated attachment exponent.\n");
  } else cat("95% confidence interval of the attachment exponent: (", object$ci[1], ",", 
       object$ci[2],")\n");
}