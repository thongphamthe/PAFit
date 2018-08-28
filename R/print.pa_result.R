# function to summarize estimation results  2015-3-11 Thong Pham
print.PA_result <- function(x,...){
  cat("\nContaining the estimated attachment function. \n");
  cat("Number of bins:", x$g,"\n");
  cat("Estimated attachment exponent:", x$alpha, "\n");
  if (x$ci[1] == "N") {
    cat("No possible confidence interval for the estimated attachment exponent.\n");
  } else cat("Attachment exponent ","\u00B1", " 2 s.d.", " : (", x$ci[1], ",", 
              x$ci[2],")\n",sep = "");
}