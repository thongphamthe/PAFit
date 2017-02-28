print.CV_Result <- function(x,...){
    cat("\nCV_Result object contains the cross validation result. \n");
    if (!is.null(x$s_optimal))
        print(paste0("Optimal r parameter is: ", x$r_optimal));
    if (!is.null(x$s_optimal))
        print(paste0("Optimal s parameter is: ", x$s_optimal));
}