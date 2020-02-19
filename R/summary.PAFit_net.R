summary.PAFit_net <- function(object,...){
  print(object)
  x <- object
  if (!is.null(x$PA)) {
    cat(paste0("\nThe PA function is: ",paste(head(x$PA), collapse = ", "),",..."))
  } else cat(paste0("\nThere is no PA function"))
  if (!is.null(x$fitness)) {
    cat(paste0("\nThe node fitnesses are: ",paste(head(x$fitness), collapse = ", "),",..."))
  } else cat(paste0("\nThere are no node fitnesses"))
  cat("\n")
}