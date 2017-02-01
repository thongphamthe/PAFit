summary.PA_result <- function(object,...){
  cat("\nPA_result object contains estimated attachment function. \n");
  cat("Number of bins: ", object$G,"\n");
  cat("Estimated attachment exponent:", object$alpha, "\n");
<<<<<<< HEAD
  if (object$ci[1] == "N") {
      cat("No possible confidence interval for the estimated attachment exponent.\n");
  } else cat("95% confidence interval of the attachment exponent: (", object$ci[1], ",", 
       object$ci[2],")\n");
=======
  cat("95% confidence interval of the attachment exponent: (", object$ci[1], ",", 
      object$ci[2],")\n");
>>>>>>> master
}