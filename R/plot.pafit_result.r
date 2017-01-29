# function to plot estimation results  2015-3-12 Thong Pham
plot.PAFit_result <-
function(x,net_stat,true_f = NULL, plot = c("A","f","true_f"), plot_bin = TRUE,
         line = FALSE, confidence = TRUE, high_deg = NULL, shade_point = 0.5, shade_interval = 0.5, 
         max_A = NULL, 
         min_A = NULL, f_min = NULL, 
         f_max = NULL, plot_true_degree = FALSE, 
         label_x = NULL, label_y = NULL,
         col_interval = "lightsteelblue",col_point = "black",...) {
  if (plot_bin == TRUE) {
      x$k <- x$center_k
      x$A <- x$theta
      x$upper_A <- x$upper_bin
      x$lower_A <- x$lower_bin
  } 
  if ("A" == plot[1]) {
      if (!is.null(high_deg))
          non_zero <- which(x$A > 10^-20 & x$k >= max(x$deg_threshold,high_deg))
      else 
          non_zero <- which(x$A > 10^-20 & x$k >= x$deg_threshold) 
      new_var_log <- x$var_logA[non_zero]*(x$A[non_zero][1]) ^ 2;
      if (!is.null(high_deg)) {
          x$A[non_zero] <- x$A[non_zero] / x$A[non_zero][1];    
          x$lower_A[non_zero] <- exp(log(x$A[non_zero]) - 2*sqrt(new_var_log));
          x$upper_A[non_zero] <- exp(log(x$A[non_zero]) + 2*sqrt(new_var_log));
      } 
     if ((!is.null(max_A)) && (!is.null(min_A)))
          limit <- c(min(min_A,x$lower_A[non_zero]),max(max_A,x$upper_A[non_zero]))
      else if (!is.null(min_A))
          limit <- c(min(min_A,x$lower_A[non_zero]),max(x$upper_A[non_zero])) 
      else if (!is.null(max_A))
          limit <- c(min(x$lower_A[non_zero]),max(max_A,x$upper_A[non_zero]))
      else 
          limit <- c(min(x$lower_A[non_zero]),max(x$upper_A[non_zero]))
    
      xlim  <- c(min(x$k[non_zero] + 1),max(x$k[non_zero] + 1))
      plot(x$k[non_zero][1] + 1,x$A[non_zero][1],xlab = ifelse(!is.null(label_x),label_x,expression(k + 1)),
           ylab = ifelse(!is.null(label_y),label_y,expression(hat(A)[k])),
           xlim = xlim, ylim = limit,log = "xy",col = col_point, ...)
      #xtick = seq(from = xlim[1], to = xlim[2],5)
      #axis(side = 1, at = xtick, labels = NULL, xlim = xlim, log = "x")
      if (TRUE == confidence) {
          order    <- order(x$k[non_zero] + 1)
          upper_f  <- x$upper_A[non_zero][order]
          lower_f  <- x$lower_A[non_zero][order]
          unique_x <- unique((x$k[non_zero] + 1)[order])
          upper_u  <- rep(0,length(unique_x))
          lower_u  <- upper_u    
          for (jjj in 1:length(unique_x)) {
            uu           <- which((x$k[non_zero] + 1)[order] == unique_x[jjj])
            #print(uu)
            upper_u[jjj] <- max(upper_f[uu])
            lower_u[jjj] <- min(lower_f[uu])
          }
          my_col <- as.vector(col2rgb(col_interval))
          my_col <- my_col / 255
          polygon(c(unique_x,rev(unique_x)),c(upper_u,rev(lower_u)), col = rgb(my_col[1],my_col[2],my_col[3],shade_interval),border = NA) 
          #arrows(x0 = x$k[non_zero] + 1, y0 = x$lower_A[non_zero], 
          #       x1 = x$k[non_zero] + 1, y1 = x$upper_A[non_zero], code = 3,angle = 90, 
          #       length = 0,
          #       col = rgb(0,0,0,shade_interval))
      }
          points(x$k[non_zero] + 1,x$A[non_zero], col = col_point,...)
      if (TRUE == line) {
          alpha <- x$alpha
          beta <-  x$loglinear_fit$coefficients[1]
          lines(x$k[non_zero],exp(beta)*(x$k[non_zero])^alpha,lwd= 2)
      }
  }
  else if ("f" == plot[1]) {
      if (FALSE == is.null(high_deg))
          non_zero <- x$lower_f > 10^-20 & net_stat$increase >= high_deg
      else
          non_zero <- x$lower_f > 10^-20 & net_stat$increase > 0
      if (length(non_zero) <= 0)
        stop("There is no data. Please decrease high_deg") 
      if (TRUE == confidence)
          lim_y = c(min(x$lower_f[non_zero]), 
                    max(x$upper_f[non_zero]))
      else lim_y = c(min(x$f[non_zero]),max(x$f[non_zero]))
      xlim <- c(min(net_stat$increase[non_zero])+1,max(net_stat$increase[non_zero]))
      plot(net_stat$increase[non_zero][1],x$f[non_zero][1],log="xy",ylab = "Estimated fitness",xlim = xlim, 
           ylim = lim_y, xlab = "Number of edges acquired",...)
      if (TRUE == confidence) {
          order   <- order(net_stat$increase[non_zero])
          upper_f <- x$upper_f[non_zero][order]
          lower_f <- x$lower_f[non_zero][order]
          unique_x <- unique(net_stat$increase[non_zero][order])
          upper_u  <- rep(0,length(unique_x))
          lower_u  <- upper_u  
        for (jjj in 1:length(unique_x)) {
          uu           <- which(net_stat$increase[non_zero][order] == unique_x[jjj])
          #print(uu)
          upper_u[jjj] <- max(upper_f[uu])
          lower_u[jjj] <- min(lower_f[uu])
        }  
          my_col <- as.vector(col2rgb(col_interval))
          my_col <- my_col / 255
          polygon(c(unique_x,rev(unique_x)),c(upper_u,rev(lower_u)), col = rgb(my_col[1],my_col[2],my_col[3],shade_interval),
                  border = NA)    
          #arrows(x0 = net_stat$increase[non_zero], y0 = x$lower_f[non_zero], x1 = net_stat$increase[non_zero], 
          #       y1 = x$upper_f[non_zero], code = 3,angle = 90, length = 0,col = rgb(0,0,0,shade_interval))
      }
      points(net_stat$increase[non_zero],x$f[non_zero],pch = 20,col = rgb(0,0,0,shade_point),...)
      abline(h = 1)
  }
  else if ("true_f" == plot[1]) {
          #names(true_f) <- net_stat$node_id 
          #true_f        <-  true_f[names(x$f)] 
          true_f1   <- length(true_f[net_stat$node_id][net_stat$f_position])*true_f[net_stat$node_id][net_stat$f_position]/
                       sum(true_f[net_stat$node_id][net_stat$f_position])
          if (FALSE == is.null(high_deg)) {
              non_zero <- x$lower_f[net_stat$f_position] > 10^-20 & true_f1 > 10^-20 & net_stat$increase[net_stat$f_position] > high_deg
          } else
              non_zero <- x$lower_f[net_stat$f_position] > 10^-20 & true_f1 > 10^-20 
          if (length(non_zero) <= 0)
             stop("There is no net_stat. Please decrease high_deg")  
          #print(non_zero)
          b        <- lm(true_f1[non_zero] ~ 0 + x$f[net_stat$f_position][non_zero])$coefficients[1]
          upper_f <-  exp(log(b*x$f[net_stat$f_position][non_zero]) + 2 * sqrt(x$var_f[net_stat$f_position][non_zero] / x$f[net_stat$f_position][non_zero] ^ 2))
          lower_f <-  exp(log(b*x$f[net_stat$f_position][non_zero]) - 2 * sqrt(x$var_f[net_stat$f_position][non_zero] / x$f[net_stat$f_position][non_zero] ^ 2))
          
          #print(lower_f)
          #print("---")
          #print(upper_f)
          if ((!is.null(f_min)) && (!is.null(f_max))) {
              if (TRUE == confidence)  
                  xlim <- c(min(c(f_min,lower_f,true_f1[non_zero])), max(c(f_max,upper_f,true_f1)))
              else xlim <- c(min(c(f_min,true_f1[non_zero])), max(c(f_max,true_f1)))
          }
          else  if (!is.null(f_max)) {
              if (TRUE == confidence)    
                  xlim <- c(min(c(lower_f,true_f1[non_zero])), max(c(f_max,upper_f,true_f1)))
              else xlim <- c(min(c(true_f1[non_zero])), max(c(f_max,true_f1)))
          }    
          else if (!is.null(f_min)) {
              if (TRUE == confidence)  
                  xlim <- c(min(c(f_min,lower_f,true_f1[non_zero])), max(c(upper_f,true_f1)))  
              else xlim <- c(min(c(f_min,true_f1[non_zero])), max(c(true_f1)))  
          }
          else {   
              if (TRUE == confidence)  
                  xlim <- c(min(c(lower_f,true_f1[non_zero])), max(c(upper_f,true_f1)))  
              else  xlim <- c(min(c(true_f1[non_zero])), max(c(true_f1)))  
          }
          #print(xlim)          
          ylim <- xlim
          plot(true_f1[non_zero][1], b*x$f[net_stat$f_position][non_zero][1], xlim= xlim, ylim = ylim,
                ylab= "Estimated fitness",xlab = "True fitness",log = "xy", pch ="",...)
          if (TRUE == confidence) {
              #x_point <- c(true_f1[non_zero],rev(true_f1[non_zero]))
              #y_point <- c(upper_f, rev(lower_f))
              order   <- order(true_f1[non_zero])
              upper_f <- upper_f[order]
              lower_f <- lower_f[order]
              unique_x <- unique(true_f1[non_zero][order])
              upper_u  <- rep(0,length(unique_x))
              lower_u  <- upper_u  
              for (jjj in 1:length(unique_x)) {
                  uu           <- which(true_f1[non_zero][order] == unique_x[jjj])
                  #print(uu)
                  upper_u[jjj] <- max(upper_f[uu])
                  lower_u[jjj] <- min(lower_f[uu])
              }
              my_col <- as.vector(col2rgb(col_interval))
              my_col <- my_col / 255
              polygon(c(unique_x,rev(unique_x)),c(upper_u,rev(lower_u)), col = rgb(my_col[1],my_col[2],my_col[3], alpha = shade_interval),border = NA)  
              #arrows(y0 = lower_f, x0 = true_f1[non_zero], y1 = upper_f, x1 = true_f1[non_zero], code = 3,angle = 90, length = 0,
              #       col = rgb(0,0,0,shade_interval))
          }
          abline(a=0,b = 1)
          if (FALSE == plot_true_degree) {
              points(true_f1[non_zero],b*x$f[net_stat$f_position][non_zero],pch = 20,col = rgb(0,0,0,shade_point),...) 
          }
          else {  
              points(b*x$f[net_stat$f_position][non_zero],true_f1[non_zero],pch ="",...)
              col <- brewer.pal(9,"Greens")
              order <- order(net_stat$increase[net_stat$f_position][non_zero])
              col_seq <- seq(min(net_stat$increase[net_stat$f_position][non_zero]),
                         max(net_stat$increase[net_stat$f_position][non_zero]),9)
              a <- sapply(net_stat$increase[net_stat$f_position][non_zero], function(x)which(col_seq >= x)[1])
              a[is.na(a)] <- 9
              text(b*x$f[net_stat$f_position][non_zero],true_f1[non_zero],net_stat$increase[net_stat$f_position][non_zero],col = col[a],...) 
          }
      }
}
