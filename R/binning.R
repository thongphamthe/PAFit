.binning <- function(deg_penmax, G = 100) {
  # deg_penmax: the penultimate degree
  # G: number of bins from 0 to deg.penmax, including deg.penmax

  start_deg <- 0
  if (G <= deg_penmax - start_deg + 1) {
    if (1 == G) {
      base <- deg_penmax - start_deg + 1
      interval_length <- deg_penmax - start_deg + 1
    }
    else {
      is.warn <- options()$warn
      options(warn = -1)
      ff <- function(x) {
        deg_penmax - start_deg + 1 - sum(floor(x^(0:(G - 
                                                    1))))
      }
      base <- uniroot(ff, interval = c(1 + 1e-15, deg_penmax - 
                                         start_deg + G + 1.1), tol = .Machine$double.eps)$root
      options(warn = is.warn)
      #print(deg_penmax - start_deg + 1)
      
      interval_length <- as.integer(floor(base^(0:(G - 1))))
      if (sum(interval_length) != deg_penmax - start_deg + 1) {
         interval_length[G] <- interval_length[G] + deg_penmax - start_deg + 1 - sum(interval_length)
      }
      #print(sum(floor(base^(0:(G - 1)))))
    }
  }
  else if ( (0 == G) || (G > deg_penmax - start_deg + 1)) {
    G <- deg_penmax - start_deg + 1
    interval_length <- rep(1, G)
    base <- 1
  }
  ######## first set #########
  bin_small   <- rep(G + start_deg - 1, deg_penmax + 1)
  b_deg_small <- c(start_deg, start_deg + cumsum(interval_length)[-G])
  e_deg_small <- b_deg_small + interval_length - 1
  if (start_deg > 0) 
    bin_small[1:start_deg] <- 0:(start_deg - 1)
  for (i in 1:G) bin_small[(b_deg_small[i]:e_deg_small[i]) + 1] <- i + start_deg - 1
  names(bin_small) <- 0:(length(bin_small) - 1) 
  if (start_deg  > 1) {
    c_bin_small        <- c(0:(start_deg - 1),sqrt(b_deg_small * e_deg_small))
  } else c_bin_small <- sqrt(b_deg_small * e_deg_small)
  
 
  return(list(bin = bin_small, center_bin = c_bin_small, 
              start = b_deg_small, end = e_deg_small, G_small = G,
              base = base, interval_length = interval_length))
}