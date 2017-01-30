print.CV_Data <- function(x,...){
  cat("\nCV_Data object contains data required in the performCV function. \n");
  cat("Number of bins: ", x$stat$G,"\n");
  cat("Ratio p:", x$p, "\n");
}