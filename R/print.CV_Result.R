print.CV_Result <- function(x,...){
    cat("\nContaining the cross validation result. \n");
    if (!is.null(x$s_optimal))
        cat(paste0("Optimal r parameter is: ", x$r_optimal, "\n"));
    if (!is.null(x$s_optimal))
        cat(paste0("Optimal s parameter is: ", x$s_optimal, "\n"));
}