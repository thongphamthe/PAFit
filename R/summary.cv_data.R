# function to summarize estimation results  2015-3-11 Thong Pham
summary.CV_Data <- function(object,...){
  cat("\nContaining the data required in cross-validation. \n");
  cat("Number of bins:", object$stat$g,"\n");
  cat("Ratio p:", object$p, "\n");
}

