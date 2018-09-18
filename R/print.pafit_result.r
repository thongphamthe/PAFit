# function to display estimation results  2015-3-11 Thong Pham
print.PAFit_result <- function(x,...) {
  cat("\nContaining the estimation results by the PAFit method. \n")
  if (x$only_PA == TRUE) {
      cat("Mode: Only the attachment function was estimated. \n")     
  } else if (x$only_f == TRUE) {
      cat("Mode: Only node fitnesses were estimated. \n")
  }
  else {
      cat("Mode: Both the attachment function and node fitness were estimated. \n")
  }
  
  #cat(" Form of the PA function:",x$mode_f,"\n");
  if (x$only_f == FALSE) {
      if (x$auto_lambda == TRUE) {
          cat("Selected r parameter:", x$ratio,"\n");  
      } else cat("Lambda used:", x$lambda,"\n");
  }
  
  if (x$only_PA == FALSE)
      cat("Selected s parameter: ",x$shape,"\n", sep = "")
  cat("Estimated attachment exponent:",x$alpha,"\n");
  if (x$ci[1] == "N") {
    cat("No possible confidence interval for the estimated attachment exponent.\n");
  } else if (x$mode_f != "Log_linear") {
        cat("Attachment exponent ","\u00B1", " 2 s.d.", ": (", x$ci[1], ",", 
           x$ci[2],")\n",sep = "");
  }
  else {
     cat("Attachment exponent ","\u00B1", " 2 s.d.", ": (", x$ci[1], ",", 
         x$ci[2],")\n",sep = "");  

  }
}