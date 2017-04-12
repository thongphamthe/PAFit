# function to summarize estimation results  2015-3-11 Thong Pham
print.PA_result <- function(x,...){
  cat("\nContaining the estimated attachment function. \n");
  cat("Number of bins: ", x$G,"\n");
  cat("Estimated attachment exponent:", x$alpha, "\n");
  if (x$ci[1] == "N") {
    cat("No possible confidence interval for the estimated attachment exponent.\n");
  } else cat("95% confidence interval of the attachment exponent: (", x$ci[1], ",", 
              x$ci[2],")\n");
}