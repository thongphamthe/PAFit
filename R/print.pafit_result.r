# function to display estimation results  2015-3-11 Thong Pham
print.PAFit_result <- function(x,...) {
  cat("\nPAFit_result object contains the estimation results by the PAFit method. \n")
  if (x$only_PA == TRUE) {
      cat("Mode: Only the PA function was estimated.")     
  } else if (x$only_f == TRUE) {
      cat("Mode: Only node fitnesses were estimated.")
  }
  else {
      cat("Mode: Both the attachment kernel and node fitness were estimated.")
  }
  cat(" Form of the PA function:",x$mode_f,"\n");
  if (x$auto_lambda == TRUE) {
    cat("Ratio (r): ", x$ratio,"\n");
  } else cat("Lambda used: ", x$lambda,"\n");
  cat("Prior of node fitness: shape: ",x$shape,"; rate: ",x$rate,"\n")
  cat("Estimated attachment exponent: ",x$alpha,"\n");
  if (x$ci[1] == "N") {
    cat("No possible confidence interval for the estimated attachment exponent.\n");
  } else if (x$mode_f != "Log_linear") {
        cat("95% confidence interval of the attachment exponent: (", x$ci[1], ",", 
           x$ci[2],")\n");
  }
  else {
     cat("Two-sigma confidence interval of the attachment exponent: (", x$ci[1], ",", 
         x$ci[2],")\n");  

  }
}