# function to plot estimation results  2015-3-12 Thong Pham
plot.PA_result <-
  function(x, net_stat,
           plot_bin = TRUE,
           high_deg = NULL,
           line = FALSE, 
           shade_point = 0.5, 
           max_A = NULL, 
           min_A = NULL, 
           label_x = NULL, 
           label_y = NULL,
           col_point = "black",...) {
    if (plot_bin == TRUE) {
      x$k <- x$center_k
      x$A <- x$theta
    } 
    if (!is.null(min_A))
        if (min_A <= 0)
            stop("min_A must be positive")
    if (!is.null(max_A))
        if (max_A <= 0)
            stop("max_A must be positive")
    if (!is.null(min_A) && (is.null(max_A)))
        if (min_A >= max_A)
            stop("max_A must be greater than min_A")
    
      if (!is.null(high_deg))
        non_zero <- which(x$A > 10^-20 & x$k >= max(high_deg,1))
      else 
        non_zero <- which(x$A > 10^-20 & x$k >= 1) 
      
      if (!is.null(high_deg)) {
        x$A[non_zero] <- x$A[non_zero] / x$A[non_zero][1];  
     
      } 
      if ((!is.null(max_A)) && (!is.null(min_A)))
          limit <- c(min(min_A, x$A[non_zero]),max(max_A,x$A[non_zero]))
      else if (!is.null(min_A))
          limit <- c(min(min_A,x$A[non_zero]),max(x$A[non_zero])) 
      else if (!is.null(max_A))
          limit <- c(min(x$A[non_zero]),max(max_A,x$A[non_zero]))
      else 
          limit <- c(min(x$A[non_zero]),max(x$A[non_zero]))
      
      xlim  <- c(min(x$k[non_zero] + 1),max(x$k[non_zero] + 1))
      plot(x$k[non_zero][1] + 1,x$A[non_zero][1],xlab = ifelse(!is.null(label_x),label_x,expression(k + 1)),
           ylab = ifelse(!is.null(label_y),label_y,expression(hat(A)[k])),
           xlim = xlim, ylim = limit,log = "xy",col = col_point, ...)
      #xtick = seq(from = xlim[1], to = xlim[2],5)
      #axis(side = 1, at = xtick, labels = NULL, xlim = xlim, log = "x")
      
      points(x$k[non_zero] + 1,x$A[non_zero], col = col_point,...)
      if (TRUE == line) {
        alpha <- x$alpha
        beta <-  x$linear_fit$coefficients[1]
        lines(x$k[non_zero],exp(beta)*(x$k[non_zero])^alpha,lwd= 2)
      }
}
