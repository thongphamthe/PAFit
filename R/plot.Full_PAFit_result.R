# function to plot estimation results  2015-3-12 Thong Pham
plot.Full_PAFit_result <-
  function(x                       ,
           net_stat                ,
           true_f         = NULL   , plot             = "A"              , plot_bin     = TRUE  ,
           line           = FALSE  , confidence       = TRUE             , high_deg_A   = 1     , 
           high_deg_f     = 5      ,   
           shade_point    = 0.5    , col_point        = "grey25"         , pch          = 16    , 
           shade_interval = 0.5    , col_interval     = "lightsteelblue" ,
           label_x        = NULL   , label_y          = NULL             ,
           max_A          = NULL   , min_A            = NULL             , f_min        = NULL  , 
           f_max          = NULL   , plot_true_degree = FALSE            , 
           ...) {
    x <- x$estimate_result
    if (plot_bin == TRUE) {
      #print("plot bin")
      x$k       <- x$center_k
      x$A       <- x$theta
      x$upper_A <- x$upper_bin
      x$lower_A <- x$lower_bin
      x$var_logA <- x$var_logbin
    } 
    #names(x$k)        <- as.character(as.integer(x$k))
    #names(x$var_logA) <- as.character(as.integer(x$k))
    #print("In plot.full_PAFit_result")
    
    dots <- function(...) {
      list(...)
    }
    additional_para <- dots(...)
    #print(additional_para)
    
    gg_color_hue = function( n ) {
      hues = seq( 15 , 375 , length = n + 1 )
      hcl( h = hues , l = 65 , c = 100 )[ 1 : n ]
    }
    cols = c( "grey25" , gg_color_hue( n = 3 ) )
    gray25 = cols[ 1 ]; red = cols[ 2 ]; green = cols[ 3 ]; blue = cols[ 4 ]
    
    if ("A" == plot[1]) {
      if (!is.null(high_deg_A))
        non_zero <- which(x$A > 10^-20 & x$k >= high_deg_A & x$k > 0)
      else 
        non_zero <- which(x$A > 10^-20 & x$k >= x$deg_threshold & x$k > 0) 
      #non_zero   <- x$k[non_zero]
      #non_zero   <- as.character(as.integer(non_zero))
      #print(x$A)
      #print(x$k)
      #print(non_zero)
      if (!is.null(high_deg_A)) {
        #print(high_deg_A)
        v                   <- which(x$k[non_zero] == high_deg_A)  
        if (length(v) > 0) {
          new_var_log         <- x$var_logA[non_zero] #* (x$A[non_zero][v[1]]) ^ 2;
          x$A[non_zero]       <- x$A[non_zero] / x$A[non_zero][v[1]];  
          #print(non_zero)
          #print(x$var_logA)
          #print(x$var_A)
          #print(x$var_logA[non_zero])
          
          #print(new_var_log)
          #print(x$var_logA[non_zero])
          #print(non_zero)
          #print(v)
          #print(x$k)
          #print(x$A)
          x$lower_A[non_zero] <- exp(log(x$A[non_zero]) - 2 * sqrt(new_var_log));
          x$upper_A[non_zero] <- exp(log(x$A[non_zero]) + 2 * sqrt(new_var_log));
          non_zero            <- x$A > 10^-20 & x$k >= high_deg_A &
            (!is.na(x$upper_A)) & (!is.infinite(x$upper_A))
        }
      } 
      if ((!is.null(max_A)) && (!is.null(min_A)))
        limit <- c(min(min_A,x$lower_A[non_zero]) , max(max_A,x$upper_A[non_zero]))
      else if (!is.null(min_A)) {
        limit <- c(min(min_A,x$lower_A[non_zero]) , max(x$upper_A[non_zero])) 
      }
      else if (!is.null(max_A)) {
        limit <- c(min(x$lower_A[non_zero]) , max(max_A,x$upper_A[non_zero]))
      }
      else { 
        limit <- c(min(x$lower_A[non_zero]) , max(x$upper_A[non_zero]))
      }
      
      if (is.null(additional_para$ylim))
        ylim <- limit  
      else ylim <- additional_para$ylim
      
      
      if (is.null(additional_para$xlim))
        xlim <- c(min(x$k[non_zero]) , max(x$k[non_zero]))  
      else xlim <- additional_para$xlim
      #print(xlim)
      
      temp   <- names(additional_para)
      ok_vec <- which(temp != "xlim" & temp != "ylim")
      #print(additional_para)
      additional_para_list <- ""
      for (i in ok_vec) {
        if (!is.character(additional_para[[i]]))  {
          additional_para_list <- paste0(additional_para_list ,temp[i]," = ",additional_para[[i]],",");
        } else  additional_para_list <- paste0(additional_para_list ,temp[i]," = '",additional_para[[i]],"',") 
      }
      final_para <- substr(additional_para_list,1,nchar(additional_para_list) - 1)
      #print(final_para)
      col_pa <- as.vector(col2rgb(col_point)) / 255
      #print(non_zero)
      eval(parse(text = paste('plot(x$k[non_zero][1], x$A[non_zero][1] , xlab = ifelse(!is.null(label_x),label_x, "Degree k"),
                              ylab = ifelse(!is.null(label_y) , label_y,"Attachment function"),
                              xlim = xlim , ylim = ylim , axes = FALSE, log = "xy" , 
                              col = rgb(col_pa[1],col_pa[2],col_pa[3], shade_point), 
                              mgp = c( 2.5 , 1 , 0 ),
                              type = "n",', final_para, ')')))
      #magaxis(grid = TRUE, frame.plot = TRUE)
      eval(parse(text = paste('magaxis(grid = FALSE, frame.plot = TRUE,',final_para,')')));
      #xtick = seq(from = xlim[1], to = xlim[2],5)
      #axis(side = 1, at = xtick, labels = NULL, xlim = xlim, log = "x")
      if (TRUE == confidence) {
        order    <- order(x$k[non_zero])
        upper_f  <- x$upper_A[non_zero][order]
        lower_f  <- x$lower_A[non_zero][order]
        unique_x <- unique((x$k[non_zero])[order])
        upper_u  <- rep(0,length(unique_x))
        lower_u  <- upper_u    
        for (jjj in 1:length(unique_x)) {
          uu           <- which((x$k[non_zero])[order] == unique_x[jjj])
          #print(uu)
          upper_u[jjj] <- max(upper_f[uu])
          lower_u[jjj] <- min(lower_f[uu])
        }
        my_col <- as.vector(col2rgb(col_interval))
        my_col <- my_col / 255
        polygon(c(unique_x,rev(unique_x)) , c(upper_u,rev(lower_u)) , col = rgb(my_col[1],my_col[2],my_col[3],shade_interval),border = NA) 
        #arrows(x0 = x$k[non_zero] + 1, y0 = x$lower_A[non_zero], 
        #       x1 = x$k[non_zero] + 1, y1 = x$upper_A[non_zero], code = 3,angle = 90, 
        #       length = 0,
        #       col = rgb(0,0,0,shade_interval))
      }
      points(x$k[non_zero], x$A[non_zero] , pch = pch, col = rgb(col_pa[1],col_pa[2],col_pa[3], shade_point),...)
      if (TRUE == line) {
        alpha <- x$alpha
        #if (!is.null(names(x$loglinear_fit)))
        #    beta <-  x$loglinear_fit$coefficients[1]
        #else {
        index_one <- which(x$A[non_zero] == 1 & x$k[non_zero] != 0)[1] 
        beta      <- -alpha* log(x$k[non_zero][index_one])
        #print(beta)
        #}
        lines(x$k[non_zero], exp(beta) * (x$k[non_zero]) ^ alpha, lwd = 2, col = blue)
        
        #lines(c(1,x$k[non_zero] + 1),c(x$PA_offset , exp(beta) * (x$k[non_zero]) ^ alpha), lwd = 2, col = green)
      }
      
      
    }
    else if ("f" == plot[1]) {
      
      non_zero <- x$lower_f > 10^-20 & net_stat$increase > 0
      if (length(non_zero) <= 0)
        stop("There is no data. Please decrease high_deg") 
      # Plot the density of node fitnesses
      
      f_non <- x$f[non_zero]
      d     <- density(f_non)
      ok_d  <- d$x > 0
      max_x <- max(d$x[ok_d])
      min_x <- min(d$x[ok_d])
      #print(min_x)
      #print(max_x)
      #layout(cbind(1,2), width = c(4,1))
      plot(d$x[ok_d] , d$y[ok_d], col  = 2 , log = "x", lwd = 0, main = "", 
           xlab = "Fitness", ylab = "Density",
           xlim = c(min_x,max_x),
           mgp = c(2.5 , 1 , 0 ), 
           axes = FALSE,
           pch = "",...)
      axis(side = 1, at=c(format(round(min_x, 2), nsmall = 2),
                          format(round(1, 2), nsmall = 2),
                          format(round(max_x, 2), nsmall = 2)),
                          tcl = 0.5, ...)
      magaxis(side = 2,grid = FALSE, frame.plot = TRUE,
              logpretty = FALSE,...)
      #x <- c(format(min(f_non),digits = 1),1,5,10,15,format(max(f_non),digits = 4))
      u    <- smooth.spline(d$x, d$y, spar = 0.01)
      #polygon(d$x, d$y, col = red_fade, border=NA)
      ok_u <- u$x > 0 & u$y > 0 
      lines(u$x[ok_u], u$y[ok_u], col = "grey50",lwd = 2.5);

    }
    else if ("true_f" == plot[1]) {
      #names(true_f) <- net_stat$node_id 
      #true_f        <-  true_f[names(x$f)] 
      true_f1   <- length(true_f[net_stat$node_id][net_stat$f_position]) * true_f[net_stat$node_id][net_stat$f_position]/
        sum(true_f[net_stat$node_id][net_stat$f_position])
      if (FALSE == is.null(high_deg_f)) {
        non_zero <- x$lower_f[net_stat$f_position] > 10^-20 & true_f1 > 10^-20 & net_stat$increase[net_stat$f_position] > high_deg_f
      } else
        non_zero <- x$lower_f[net_stat$f_position] > 10^-20 & true_f1 > 10^-20 
      if (length(non_zero) <= 0)
        stop("There is no net_stat. Please decrease high_deg")  
      #print(non_zero)
      b        <- lm(true_f1[non_zero] ~ 0 + x$f[net_stat$f_position][non_zero])$coefficients[1]
      upper_f  <-  exp(log(b*x$f[net_stat$f_position][non_zero]) + 2 * sqrt(x$var_f[net_stat$f_position][non_zero] / x$f[net_stat$f_position][non_zero] ^ 2))
      lower_f  <-  exp(log(b*x$f[net_stat$f_position][non_zero]) - 2 * sqrt(x$var_f[net_stat$f_position][non_zero] / x$f[net_stat$f_position][non_zero] ^ 2))
      
      #print(lower_f)
      #print("---")
      #print(upper_f)
      
      if (is.null(additional_para$xlim)) {
        if ((!is.null(f_min)) && (!is.null(f_max))) {
          if (TRUE == confidence)  
            xlim <- c(min(c(f_min,lower_f,true_f1[non_zero])) , max(c(f_max,upper_f,true_f1)))
          else xlim <- c(min(c(f_min,true_f1[non_zero])) , max(c(f_max,true_f1)))
        }
        else  if (!is.null(f_max)) {
          if (TRUE == confidence)    
            xlim <- c(min(c(lower_f,true_f1[non_zero])) , max(c(f_max,upper_f,true_f1)))
          else xlim <- c(min(c(true_f1[non_zero])) , max(c(f_max,true_f1)))
        }    
        else if (!is.null(f_min)) {
          if (TRUE == confidence)  
            xlim <- c(min(c(f_min,lower_f,true_f1[non_zero])) , max(c(upper_f,true_f1)))  
          else xlim <- c(min(c(f_min,true_f1[non_zero])) , max(c(true_f1)))  
        }
        else {   
          if (TRUE == confidence)  
            xlim <- c(min(c(lower_f,true_f1[non_zero])) , max(c(upper_f,true_f1)))  
          else  xlim <- c(min(c(true_f1[non_zero])) , max(c(true_f1)))  
        }
      }
      else {
        xlim <- additional_para$xlim  
      }
      #print(xlim)  
      if (is.null(additional_para$ylim)) 
        ylim <- xlim
      else ylim <- additional_para$ylim
      temp   <- names(additional_para)
      ok_vec <- which(temp != "xlim" & temp != "ylim")
      #print(additional_para)
      additional_para_list <- ""
      for (i in ok_vec) {
        if (!is.character(additional_para[[i]]))  {
          additional_para_list <- paste0(additional_para_list ,temp[i]," = ",additional_para[[i]],",");
        } else  additional_para_list <- paste0(additional_para_list ,temp[i]," = '",additional_para[[i]],"',") 
      }
      final_para <- substr(additional_para_list,1,nchar(additional_para_list) - 1)
      
      #print(final_para)
      
      eval(parse(text = paste('plot(true_f1[non_zero][1], b * x$f[net_stat$f_position][non_zero][1] , xlim= xlim , 
                              ylim = ylim,
                              ylab = "Estimated fitness" , xlab = "True fitness" , log = "xy" , type = "n", axes = FALSE,
                              ,
                              mgp = c( 2.5 , 1 , 0 ) ,', 
                              final_para, ')')))
      #magaxis(grid = TRUE, frame.plot = TRUE)
      eval(parse(text = paste('magaxis(grid = FALSE, frame.plot = TRUE,',final_para,')')));
      
      
      if (TRUE == confidence) {
        #x_point <- c(true_f1[non_zero],rev(true_f1[non_zero]))
        #y_point <- c(upper_f, rev(lower_f))
        order    <- order(true_f1[non_zero])
        upper_f  <- upper_f[order]
        lower_f  <- lower_f[order]
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
        polygon(c(unique_x , rev(unique_x)) , c(upper_u , rev(lower_u)) , col = rgb(my_col[1] , my_col[2] , my_col[3] , alpha = shade_interval) , border = NA)  
        #arrows(y0 = lower_f, x0 = true_f1[non_zero], y1 = upper_f, x1 = true_f1[non_zero], code = 3,angle = 90, length = 0,
        #       col = rgb(0,0,0,shade_interval))
      }
      abline(a=0,b = 1)
      if (FALSE == plot_true_degree) {
        col_fit <- as.vector(col2rgb(col_point)) / 255
        points(true_f1[non_zero],b * x$f[net_stat$f_position][non_zero] , pch = pch , col = rgb(col_fit[1],col_fit[2],
                                                                                                col_fit[3],
                                                                                                shade_point),...) 
      }
      else {  
        points(b * x$f[net_stat$f_position][non_zero],true_f1[non_zero] , pch = "" , ...)
        col         <- brewer.pal(9 , "Greens")
        order       <- order(net_stat$increase[net_stat$f_position][non_zero])
        col_seq     <- seq(min(net_stat$increase[net_stat$f_position][non_zero]) ,
                           max(net_stat$increase[net_stat$f_position][non_zero]) , 9)
        a           <- sapply(net_stat$increase[net_stat$f_position][non_zero] , function(x)which(col_seq >= x)[1])
        a[is.na(a)] <- 9
        text(b * x$f[net_stat$f_position][non_zero] , true_f1[non_zero] , net_stat$increase[net_stat$f_position][non_zero] , col = col[a] , ...) 
      }
    }
  }
