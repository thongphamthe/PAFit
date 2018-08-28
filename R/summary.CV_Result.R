summary.CV_Result <- function(object,...){
    cat("\nContaining the cross validation result. \n");
    if (!is.null(object$r_optimal))
        cat(paste0("Selected r parameter is:",object$r_optimal, "\n"));
    if (!is.null(object$s_optimal))
        cat(paste0("Selected s parameter is:",object$s_optimal, "\n"));
}