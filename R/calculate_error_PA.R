

calculate_error_PA <- function(k , A , start_deg = 0 , mode = 1 , alpha = 1 , beta = 1 , sat_at = 100) {
  #nonzero <- 0:(max(max.deg-1,k))
  
  nonzero <- k > start_deg
  A       <- A[nonzero]
  k       <- k[nonzero]
  nonzero <- k
  if (mode == 1) {          #BA model
      true <- pmax(nonzero,1) ^ alpha
  } else if (mode == 2) {     #saturation model
      true <- (pmin(nonzero,sat_at))^alpha
  }else if (mode == 3)     #log linearity
      true <- alpha*log(pmax(nonzero,1))^beta + 1
  
  temp <- lm(true ~ 0 + A, weights = 1 / true ^ 2) 
  return(sum(1 / length(true) * temp$residual ^ 2 * temp$weights))
}

