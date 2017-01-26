# function to display estimation results  2015-3-11 Thong Pham
print.PAFit_result <- function(x,...) {
  cat("\nPAFit_result object contains the estimation results by the PAFit method. \n")
  if (x$only_PA == TRUE) {
      cat("Mode: Only the attachment kernel was estimated.")
  }
  else if (x$only_f == TRUE) {
    cat("Mode: Only node fitnesses were estimated.")
  }
  else {
    cat("Mode: Both the attachment kernel and node fitness were estimated.")
  }
}