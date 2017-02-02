# function to generate generalized Erdos-Renyi network 
Generate_BB <-
  function(N, 
           num_seed       = 2      , 
           multiple_node  = 1      , 
           m              = 1      ,
           mode_f         = "gamma", 
           rate           = 10     , 
           shape          = 10     , 
           meanlog        = 0      , 
           sdlog          = 1      ,
           scale_pareto   = 2      ,
           shape_pareto   = 2      
  ){
    # N: number of nodes
    # Number of time-step: (N  - num_seed) / multiple_node
    if (num_seed >= N)
        stop("num_seed too large")
    if (num_seed < 2)
        stop("num_seed too small")
    if (multiple_node > N - num_seed)
        stop("Multiple node and/or num_seed are too large")   
    if ((m <= 0) || (num_seed <= 0) || (multiple_node <= 0) || (rate < 0) || (shape <= 0) || (scale_pareto <= 0) || (shape_pareto <= 0))
        stop("The parameters must be positive") 
    if (sdlog <= 0)
        stop("sdlog must be positive")   
    
    return(GenerateNet(N = N , num_seed = num_seed , multiple_node = multiple_node , m = m , alpha = 1, mode_f = mode_f, 
                       rate = rate, shape = shape,meanlog = meanlog, sdlog = sdlog, scale_pareto = scale_pareto, shape_pareto = shape_pareto))
  }
