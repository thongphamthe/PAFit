
## Assume that A is normalized at A_0 or A_1, i.e., no interception on log scale
.regress_alpha <- function(k,A,sd_log_A = NULL) {
  # upper, lower: 2 sd in log scale
  # regress the form (max(k,1))^alpha
  k_plus   <- k 
  limit    <- 10^-12
  #print(k)
  #print(A)
  ok_index <- A > limit & k_plus > limit
  x        <- log(k_plus[ok_index])
  y        <- log(A[ok_index])

  if (!is.null(sd_log_A) && (sum(is.na(sd_log_A)) < length(sd_log_A))) {
      #print("In here")
    A_ok_var_zero <- sd_log_A[ok_index] < limit
    A_ok_var_ok   <- sd_log_A[ok_index] > limit
    if (sum(A_ok_var_ok) > sum(A_ok_var_zero)) {
      var_log <- sd_log_A[ok_index][A_ok_var_ok]^2
      more_ok <- var_log > 0 & !is.nan(var_log)
      if (sum(more_ok) > 1) {
          x <- x[A_ok_var_ok][more_ok]
          y <- y[A_ok_var_ok][more_ok]
          var_log <- var_log[more_ok]
          result  <- lm(y ~ x,weights = 1/var_log)
      } else {
        result  <- lm(y ~ x)
      }
    } else {result  <- lm(y ~ x)}
  } else {
    result  <- lm(y ~ x)
  }
  names(result$coefficients) <- c("Factor","Attachment exponent")
  res                        <- df.residual(result)
  if (res > 0) {
    ci            <- confint(result,"Attachment exponent")
  } else ci <- "NA"
  alpha           <- result$coefficients[2]
  return(list(alpha = alpha, ci = ci))
}
