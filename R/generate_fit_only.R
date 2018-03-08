# function to generate generalized Erdos-Renyi network 
generate_fit_only <-
  function(N              = 1000   , 
           num_seed       = 2      , 
           multiple_node  = 1      , 
           m              = 1      ,
           mode_f         = "gamma", 
           s              = 10    
  ){
    # N: number of nodes
    # Number of time-step: (N  - num_seed) / multiple_node
    if (num_seed >= N)
        stop("num_seed too large")
    if (num_seed < 2)
        stop("num_seed too small")
    if (multiple_node > N - num_seed)
        stop("Multiple node and/or num_seed are too large")   
    if ((m <= 0) || (num_seed <= 0) || (multiple_node <= 0) || (s < 0))
        stop("The parameters must be positive") 
    
    return(generate_net(N = N , num_seed = num_seed , multiple_node = multiple_node , m = m , alpha = 0, mode_f = mode_f, 
                       s = s))
  }
