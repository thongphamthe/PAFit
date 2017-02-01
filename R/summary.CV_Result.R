summary.CV_Result <- function(object,...){
    cat("\nCV_Result object contains the cross validation result. \n");
    print(paste0("Optimal r parameter is: ",object$r_optimal));
    print(paste0("Optimal s parameter is: ",object$s_optimal));
}