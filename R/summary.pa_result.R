summary.PA_result <- function(object,...){
  cat("\nPA_result object contains estimated attachment function. \n");
  cat("Number of bins: ", object$G,"\n");
  cat("Estimated attachment exponent:", object$alpha, "\n");
}