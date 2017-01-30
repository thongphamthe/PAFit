# function to summarize estimation results  2015-3-11 Thong Pham
summary.CV_Data <- function(object,...){
  cat("\nCV_Data object contains data required in the performCV function. \n");
  cat("Number of bins: ", object$stat$G,"\n");
  cat("Ratio p:", object$p, "\n");
}

