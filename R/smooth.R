.smooth <- function(series, variance, window = round(length(series[series > -1])/2),threshold = 0.1) {
  
  # some sanity check at the start
  if (sum(series > -1) == 0) {
    print("Zero series")
    return(series)
  }
  #check minus 1 in the middle of the series
  minus_index <- which(series == -1)
  if (length(minus_index) > 0) {
      if (sum(diff(minus_index)) != (length(minus_index) - 1)) {
          stop("Middle minus 1.");  
      }
  }
  series_old <- series
  
  variance   <- variance[series > -1]
  series     <- series[series > -1]
        
  n          <- length(series)
  result     <- series
  var        <- variance
  for (j in 1:n) {
      start_index <- j - window
      end_index   <- j + window
      if (start_index <= 0) start_index <- 1;
      if (end_index > n) end_index <- n;
      displace_vector <- start_index:end_index - j
      zero_var <- which(variance[start_index:end_index] == 0)
      if (length(zero_var) > 0) {
        stop("Zero variance?")
      }
      h_value  <- seq(0.1,100,length.out = 1000)
      final_sd <- rep(0,length(h_value))
      flag     <- FALSE
      for (ii in 1:length(h_value)) {  
          weight            <- dnorm(displace_vector, mean = 0, sd = h_value[ii], log = FALSE)
          
          #min_expo          <- mean(weight[-zero_var])
          #weight[-zero_var] <- weight[-zero_var] - min_expo
          if (length(zero_var) > 0) {
              weight[-zero_var] <- weight[-zero_var]/sum(weight[-zero_var])
          } else {
              weight <- weight /sum(weight)  
          }
          if (length(zero_var) > 0)
               weight[zero_var]  <- 0
          if (length(weight) != length(series[start_index:end_index])) {
              print(paste0("Length weight: ",length(weight),"; Length series:",length(series[start_index:end_index])))
              print(start_index)
              print(end_index)
              stop("Error!")
          }
          result[j] <- sum(weight * series[start_index:end_index]) 
          var[j]    <- sum(weight^2 * variance[start_index:end_index])     
      if (sqrt(var[j]) / result[j] <= threshold) {
          flag <- TRUE
          break;  
      } else { 
        final_sd[ii] <- sqrt(var[j]) / result[j]  
      } 
      }
      if (flag == FALSE) {
         min_index <- which.min(final_sd)
         h_ok      <- h_value[min_index]
         weight           <- dnorm(displace_vector, mean = 0, sd = h_ok , log = FALSE)
         if (length(zero_var) > 0)
             weight[zero_var] <- 0
         weight <- weight/sum(weight)
         result[j] <- sum(weight * series[start_index:end_index]) 
        #print(paste0("Not found suitable h at time ",j," for the threshold ",threshold));
      }
  }
  if (length(minus_index) > 0)
      series_old[-minus_index] <- result
  return(series_old)
}