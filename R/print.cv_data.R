print.CV_Data <- function(x,...){
  cat("\nContaining the data required in cross validation. \n");
  cat("Number of bins:", x$stat$g,"\n");
  cat("Ratio p:", x$p, "\n");
}