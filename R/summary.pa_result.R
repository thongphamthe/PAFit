summary.PA_result <- function(object,...){
  cat("\nContaining the estimated attachment function. \n");
  cat("Number of bins:", object$g,"\n");
  cat("Estimated attachment exponent:", object$alpha, "\n");
  if (object$ci[1] == "N") {
      cat("No possible confidence interval for the estimated attachment exponent.\n");
  } else cat("Attachment exponent ","\u00B1", " 2 s.d.", ": (", object$ci[1], ",", 
       object$ci[2],")\n",sep = "");
}