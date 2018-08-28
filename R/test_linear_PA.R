### this function compares a set of limiting degree distribution to find the best fit to the data
### the set of limiting degree distribution includes 5 distributions:
### Yule distribution
### Waring distribution
### Negative binomial distribution
### Geometric distribution
### Poisson distribution


test_linear_PA <- function(degree_vector) {
  yule_est        <- function(degree_dist, k_min = 0) {
      n             <- sum(degree_dist)  
     
      degree        <- as.integer(labels(degree_dist)[[1]])
      upper_tail    <- degree_dist[degree > k_min]
      upper_deg     <- degree[degree > k_min]
      lower_tail    <- degree_dist[degree <= k_min]
      lower_deg     <- degree[degree <= k_min]
      dim           <- length(lower_deg)
      
      log_likelihood <- function(p = 3) {
        total_ll      <- 0  
        temp          <- 0
        if (length(lower_deg) > 0)
            for (d in 1:length(lower_deg)) {
                total_ll <- total_ll + lower_tail[d]* log(lower_tail[d] / n)
                temp     <- temp + 1 - lower_tail[d] / n
            }
        total_ll <- total_ll + (n - sum(lower_tail)) * log(temp)
        if (length(upper_deg) > 0)
            for (d in 1:length(upper_deg))  
                total_ll <- total_ll + upper_tail[d] * (log(p - 1) + lgamma(upper_deg[d]) + lgamma(p) - lgamma(upper_deg[d] + p))
        return(total_ll)
      }
      if (length(upper_deg) > 0) {
          optim_result <- optim(par = 3,fn = log_likelihood, method = "L-BFGS",
                                control = list(fnscale = -1), lower = 2 + 10^-10, upper = Inf)
          result <- list(par = optim_result$par, val = optim_result$value, AIC = -2*optim_result$value + 2*dim,
                     BIC = -2*optim_result$value + log(n) * dim, k_min = k_min, dim = dim)  
      } else {
        result <- list(par = c(NULL), val = log_likelihood(0), AIC = -2*log_likelihood(0) + 2*dim,
                       BIC = -2*log_likelihood(0) + log(n) * dim, k_min = k_min, dim = dim)    
      }
}
  
  waring_est      <- function(degree_dist, k_min = 0) {
      n             <- sum(degree_dist)  
      degree        <- as.integer(labels(degree_dist)[[1]])
      upper_tail    <- degree_dist[degree > k_min]
      upper_deg     <- degree[degree > k_min]
      lower_tail    <- degree_dist[degree <= k_min]
      lower_deg     <- degree[degree <= k_min]
      dim           <- length(lower_deg) + 1
     
      log_likelihood <- function(par = c(3,1)) {
          p             <- par[1]
          alpha         <- par[2]
          total_ll      <- 0  
          temp          <- 0
          #print("Reached here")
          if (length(lower_deg) > 0) {
              for (d in 1:length(lower_deg)) {
                  total_ll <- total_ll + lower_tail[d]* log(lower_tail[d] / n)
                  temp     <- temp + 1 - lower_tail[d] / n
              }
              total_ll <- total_ll + (n - sum(lower_tail)) * log(temp)
          }
          if (length(upper_deg) > 0) {
              for (d in 1:length(upper_deg))  
                  total_ll <- total_ll + upper_tail[d] * (log(p - 1) + lgamma(p + alpha) - lgamma(alpha + 1) + lgamma(upper_deg[d] + alpha) 
                                                          - lgamma(upper_deg[d] + p + alpha))
          }
          #print("Reached here here")
          return(total_ll)
      }
    if (length(upper_deg) > 0) {
        optim_result <- optim(par = c(3,1),fn = log_likelihood, method = "L-BFGS",control = list(fnscale = -1),
                              lower = c(2 + 10^-10,-1 + 10^-10), upper = c(Inf,Inf))
        result <- list(par = optim_result$par, val = optim_result$value, AIC = -2*optim_result$value + 2*dim,
                   BIC = -2*optim_result$value + log(n) * dim, k_min = k_min, dim = dim)  
    } else {
        result <- list(par = c(NULL,NULL), val = log_likelihood(), AIC = -2*log_likelihood() + 2*dim,
                       BIC = -2*log_likelihood() + log(n) * dim, k_min = k_min, dim = dim)    
    }
  }
  
  geom_est      <- function(degree_dist, k_min = 0) {
    n             <- sum(degree_dist)  
    degree        <- as.integer(labels(degree_dist)[[1]])
    upper_tail    <- degree_dist[degree > k_min]
    upper_deg     <- degree[degree > k_min]
    lower_tail    <- degree_dist[degree <= k_min]
    lower_deg     <- degree[degree <= k_min]
    dim           <- length(lower_deg)
    
    log_likelihood <- function(p = 0.5) {
      total_ll      <- 0  
      temp          <- 0
      if (length(lower_deg) > 0)
          for (d in 1:length(lower_deg)) {
              total_ll <- total_ll + lower_tail[d]* log(lower_tail[d] / n)
              temp     <- temp + 1 - lower_tail[d] / n
          }
      total_ll <- total_ll + (n - sum(lower_tail)) * log(temp)
      if (length(upper_deg) > 0)
          for (d in 1:length(upper_deg))  
              total_ll <- total_ll + upper_tail[d] * (log(p) + (upper_tail[d] - 1) * log(1 - p))
      return(total_ll)
    }
    if (length(upper_deg) > 0) {
        optim_result <- optim(0.5,log_likelihood, method = "L-BFGS",control = list(fnscale = -1),
                              lower = 10^-10, upper = 1 - 10^-10)
        result <- list(par = optim_result$par, val = optim_result$value, AIC = -2*optim_result$value + 2*dim,
                   BIC = -2*optim_result$value + log(n) * dim, k_min = k_min, dim = dim)  
    } else {
        result <- list(par = c(NULL), val = log_likelihood(), AIC = -2*log_likelihood() + 2*dim,
                     BIC = -2*log_likelihood() + log(n) * dim, k_min = k_min, dim = dim)  
    }
  }
  
  poisson_est   <- function(degree_dist, k_min = 0) {
    n             <- sum(degree_dist)  
   
    degree        <- as.integer(labels(degree_dist)[[1]])
    upper_tail    <- degree_dist[degree > k_min]
    upper_deg     <- degree[degree > k_min]
    lower_tail    <- degree_dist[degree <= k_min]
    lower_deg     <- degree[degree <= k_min]
    dim           <- length(lower_deg)
    
    log_likelihood <- function(p = 0.5) {
      total_ll      <- 0  
      temp          <- 0
      if (length(lower_deg) > 0)
          for (d in 1:length(lower_deg)) {
              total_ll <- total_ll + lower_tail[d]* log(lower_tail[d] / n)
              temp     <- temp + 1 - lower_tail[d] / n
          }
      total_ll <- total_ll + (n - sum(lower_tail)) * log(temp)
      if (length(upper_deg) > 0)
          for (d in 1:length(upper_deg))  
              total_ll <- total_ll + upper_tail[d] * (-p + upper_deg[d] * log(p) - lfactorial(upper_deg[d]))
      return(total_ll)
    }
    if (length(upper_deg) > 0) {
        optim_result <- optim(1,log_likelihood, method = "L-BFGS",control = list(fnscale = -1),
                              lower = 0 + 10^-10, upper = Inf)
        result <- list(par = optim_result$par, val = optim_result$value, AIC = -2*optim_result$value + 2*dim,
                       BIC = -2*optim_result$value + log(n) * dim, k_min = k_min, dim = dim)  
    } else {
        result <- list(par = c(NULL), val = log_likelihood(), AIC = -2*log_likelihood() + 2*dim,
                       BIC = -2*log_likelihood() + log(n) * dim, k_min = k_min, dim = dim)  
    }
  }
  
  nb_est         <- function(degree_dist, k_min = 0) {
    n             <- sum(degree_dist)  
    
    degree        <- as.integer(labels(degree_dist)[[1]])
    upper_tail    <- degree_dist[degree > k_min]
    upper_deg     <- degree[degree > k_min]
    lower_tail    <- degree_dist[degree <= k_min]
    lower_deg     <- degree[degree <= k_min]
    dim           <- length(lower_deg) + 1  # 2 parameter for the nb
    
    log_likelihood <- function(par = c(0.5,0.5)) {
      p             <- par[1]
      r             <- par[2]
      total_ll      <- 0  
      temp          <- 0
      if (length(lower_deg) > 0)
          for (d in 1:length(lower_deg)) {
              total_ll <- total_ll + lower_tail[d]* log(lower_tail[d] / n)
              temp     <- temp + 1 - lower_tail[d] / n
          }
      total_ll <- total_ll + (n - sum(lower_tail)) * log(temp)
      if (length(upper_deg) > 0)
          for (d in 1:length(upper_deg))  
              total_ll <- total_ll + upper_tail[d] * (lgamma(upper_deg[d] + r) - lfactorial((upper_deg[d])) - lgamma(r) +
                                                upper_deg[d] * log(p) + r * log(1 - p))
        
      return(total_ll)
    }
    if (length(upper_deg) > 0) {
        optim_result <- optim(c(0.5,1),log_likelihood, method = "L-BFGS",control = list(fnscale = -1),
                              lower = c(10^-10,10^-10),upper = c(1-10^-10,Inf))
        result <- list(par = optim_result$par, val = optim_result$value, AIC = -2*optim_result$value + 2*dim,
                   BIC = -2*optim_result$value + log(n) * dim, k_min = k_min, dim = dim)  
    } else {
      result <- list(par = c(NULL,NULL), val = log_likelihood(), AIC = -2*log_likelihood() + 2*dim,
                     BIC = -2*log_likelihood() + log(n) * dim, k_min = k_min, dim = dim)    
    }
  }
  

  
  ###############################################
  ########## main body #############
  degree_dist   <- table(degree_vector)
  degree        <- as.integer(labels(degree_dist)[[1]])
  K_max         <- max(degree_vector)
  k_min         <- degree
  full_result   <- vector(mode = "list", length = 5) # five models
  names(full_result) <- c("yule","waring","nb","geom","pois")
  for (l in 1:length(full_result))
      full_result[[l]] <- vector(mode = "list", length = length(degree))
  for (i in 1:length(k_min)) {
      yule_result      <- yule_est(degree_dist, k_min[i])
      waring_result    <- waring_est(degree_dist, k_min[i])
      nb_result        <- nb_est(degree_dist, k_min[i])
      geom_result      <- geom_est(degree_dist, k_min[i])
      poisson_result   <- poisson_est(degree_dist,k_min[i])
      full_result[[1]][[i]]      <- yule_result 
      names(full_result[[1]])[i] <- k_min[i]
      full_result[[2]][[i]]      <- waring_result 
      names(full_result[[2]])[i] <- k_min[i]
      full_result[[3]][[i]]      <- nb_result 
      names(full_result[[3]])[i] <- k_min[i]
      full_result[[4]][[i]]      <- geom_result 
      names(full_result[[4]])[i] <- k_min[i]
      full_result[[5]][[i]]      <- poisson_result 
      names(full_result[[5]])[i] <- k_min[i]
  } 
  class(full_result) <- "Linear_PA_test_result"
  return(full_result)
} 