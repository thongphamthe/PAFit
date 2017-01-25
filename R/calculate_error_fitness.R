#average error in estimating f
calculate_error_fitness <- function(true, estimate) {
  nonzero <- true > 0  
  true    <- true[nonzero]
  f       <- f[nonzero]
  n       <- length(f)
  temp    <- lm(true ~ 0 + f, weights = 1/true^2)
  return(sum(1/n*temp$residual^2*temp$weights))
}