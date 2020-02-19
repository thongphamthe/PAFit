# function to plot estimation results  2015-3-12 Thong Pham
plot.PA_result <-
  function(x, net_stat,
           plot_bin    = TRUE   ,
           high_deg    = 1      ,
           line        = FALSE  , 
           col_point   = "black",
           shade_point = 0.5    , 
           pch         = 16     ,
           max_A       = NULL   , 
           min_A       = NULL   , 
           label_x     = NULL   , 
           label_y     = NULL   ,
           ...) {
    if (plot_bin == TRUE) {
      x$k <- x$center_k
      x$A <- x$theta
      x$var_logA <- x$var_logbin
    } 
    
    dots <- function(...) {
      list(...)
    }
    additional_para <- dots(...)
    
    
    gg_color_hue = function( n ) {
      hues = seq( 15 , 375 , length = n + 1 )
      hcl( h = hues , l = 65 , c = 100 )[ 1 : n ]
    }
    cols = c( "grey25" , gg_color_hue( n = 3 ) )
    gray25 = cols[ 1 ]; red = cols[ 2 ]; green = cols[ 3 ]; blue = cols[ 4 ]
    
    
    
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
        v                   <- which(x$k[non_zero] == high_deg)  
        if (length(v) > 0) {
            x$A[non_zero]       <- x$A[non_zero] / x$A[non_zero][v[1]];  
        }
      } 
      if ((!is.null(max_A)) && (!is.null(min_A)))
          limit <- c(min(min_A, x$A[non_zero]) , max(max_A,x$A[non_zero]))
      else if (!is.null(min_A))
          limit <- c(min(min_A,x$A[non_zero]) , max(x$A[non_zero])) 
      else if (!is.null(max_A))
          limit <- c(min(x$A[non_zero]) , max(max_A,x$A[non_zero]))
      else 
          limit <- c(min(x$A[non_zero]) , max(x$A[non_zero]))
      
      if (is.null(additional_para$ylim))
          ylim <- limit  
      else ylim <- additional_para$ylim
      
      if (is.null(additional_para$xlim))
          xlim <- c(min(x$k[non_zero]) , max(x$k[non_zero]))  
      else xlim <- additional_para$xlim
      
    col_pa <- as.vector(col2rgb(col_point)) / 255
    #print("In plot.pa")
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
    
    eval(parse(text = paste('plot(x$k[non_zero][1], x$A[non_zero][1] , xlab = ifelse(!is.null(label_x),label_x, "Degree k"),
           ylab = ifelse(!is.null(label_y) , label_y,"Attachment function"), pch =pch,  mgp = c( 2.5 , 1 , 0 ),
           xlim = xlim , ylim = ylim , log = "xy" , axes = FALSE, col = col_point, type = "n",', final_para, ')')))
    eval(parse(text = paste('magaxis(grid = FALSE, frame.plot = TRUE, ',final_para,')')));
      #xtick = seq(from = xlim[1], to = xlim[2],5)
      #axis(side = 1, at = xtick, labels = NULL, xlim = xlim, log = "x")
      
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
          lines(x$k[non_zero] , exp(beta) * (x$k[non_zero]) ^ alpha , lwd= 2, col = blue)
      }
}
