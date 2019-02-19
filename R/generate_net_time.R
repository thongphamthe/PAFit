# function to generate simulated network  with PA and node fitness changes over time
.generate_net_time <-
  function(num_seed           = 100   , # number of nodes at the initial stage
           num_period         = 50    , # number of periods
           period_length      = 10    , # number of time-steps inside each period
           multiple_node      = 1     , # number of new nodes appears at each time-steps
           m                  = 1     , # number of new edges per new nodes at each time-step
           alpha_start        = 0.25  , # starting value of alpha at the start of the period
           alpha_end          = 0.75  , # ending value of alpha at the end of the period
           s                  = 10    , # inverse variance of gamma where new node fitnesses will be drawn from
           h                  = 0.1     # s.d. of the gaussian random walk at each time-step
  ) {
 
    # all fitness first value
    fit_all_first_value <- rgamma(n = num_seed + multiple_node * period_length * num_period,
                                  shape = s, rate = s)
    names(fit_all_first_value) <- 1:(num_seed + multiple_node * period_length * num_period)
   
    n_current          <- num_seed       # current number of nodes
    
    # indicates the birth period of each node
    birth_time        <- rep(0,num_seed)
    names(birth_time) <- as.numeric(1:num_seed)
    
    # The easiest way is to model each as randomly move up/down on log-scale with some factors
    # the multiplicative factor is sample from a log-normal distribution with mean 1 and sd h
    sd        <- rep(h,length = num_seed)
    names(sd) <- as.numeric(1:num_seed)

    
    fit_lifetime       <-  vector(mode = "list",length = num_period) # the lifetime fitness
    # we have to normalize fit_initial, taking into account on the node fitness that appear in that period
    #mul_factor         <-  sum(fit_all_first_value[num_seed + multiple_node * period_length]) / (num_seed + multiple_node * period_length)
    fit_initial_full   <- fit_all_first_value[1:(num_seed + multiple_node * period_length)]
    mul_factor         <- mean(fit_initial_full)
    fit_lifetime[[1]]  <- fit_initial_full / mul_factor
    
    
    alpha_lifetime     <- seq(alpha_start,alpha_end,length.out = num_period)
    
    #print(paste0("True alpha: ",alpha_lifetime))
    
    # number of rows in the edge matrix
    num_of_edge <- (num_seed - 1) + num_period * period_length * multiple_node *  m
    edgelist    <- matrix(0,nrow = num_of_edge,ncol = 3)
    
    for (i in 2:num_seed) {
      edgelist[i - 1,] <- c(i,i - 1,0);   
    }
    row_current <- num_seed - 1   # the current row index in the edgelist matrix
    
    degree                 <- rep(0,num_seed)
    names(degree)          <- as.numeric(1:num_seed) 
    u                      <- table(edgelist[1:row_current,2])
    degree[as.character(as.numeric(labels(u)[[1]]))] <- degree[as.character(as.numeric(labels(u)[[1]]))] + u 
    #print(degree)
    
    P                      <- degree
    P[P == 0]              <- 1
    
    current_time_step      <- 0
    n_old                  <- n_current
    
    for (p in 1:num_period) {
        #print(p)
        if ( p > 1) { # then changes fitness at this period p
         #print("P > 1")
          fit_all_first_value <- fit_all_first_value * rlnorm(n = length(fit_all_first_value), meanlog = 0, sdlog = h)
          new_fit_full        <- fit_all_first_value[1:(num_seed + p * multiple_node * period_length)]
          mul_factor          <- mean(new_fit_full)
          fit_lifetime[[p]]   <- new_fit_full / mul_factor
        }
        temp_fit      <- unlist(fit_lifetime[[p]])
        
        for (t in 1:period_length){
            fitness           <- temp_fit[1:n_current]
            n_old             <- n_current
            P                 <- degree
            P[P == 0]         <- 1
            current_time_step <- current_time_step + 1
            PA            <- pmax(P,1)^alpha_lifetime[p]
            PA[PA == 0]   <- 1
            node.weights  <- PA*fitness/sum(PA*fitness)
            # adding new nodes
            for (i in 1:multiple_node){
              n_current <- n_current + 1
              nodes     <- sort(sample(1:(n_old),size = m,prob = node.weights, replace = TRUE)) 
              
              ##########################################
              old_names      <- names(degree)
              degree         <- c(degree,0)
              names(degree)  <- c(old_names,as.character(as.numeric(n_current)))
              if (0 != length(nodes)) {
                temp    <- table(nodes)
                #print(paste0("Nodes:",nodes))
                node_name <- as.numeric(labels(temp)[[1]])
                degree[as.character(node_name)] <- degree[as.character(node_name)] + temp
                # for(i in 1:length(temp)) { 
                #   num_edge            <- as.numeric(temp[i]) 
                #   node_name           <- as.numeric(labels(temp[i]))
                #   print(paste0("Names:",node_name))
                #   degree[as.character(node_name)]   <- degree[as.character(node_name)] + num_edge # Update degrees.
                # }
              }
              if (0 != length(nodes)) {
                if (row_current + length(nodes) > dim(edgelist)[1]) {
                  print("Over bound")
                  edgelist           <- rbind(edgelist , matrix(0,nrow = row_current + length(nodes) -  dim(edgelist)[1], ncol = 3))
                }
              }
              
              if (0 != length(nodes)) {
                edgelist[(row_current + 1):(row_current + length(nodes)),] <- cbind(rep(n_current, length(nodes)),nodes, 
                                                                                    rep(current_time_step,length(nodes)))
              }
              row_current     <- row_current + length(nodes)
            }
            n_old             <- n_current
            P                 <- degree
            P[P == 0]         <- 1
        }  
        fit_lifetime[[p]] <- fitness
    }
    if (row_current != dim(edgelist)[1]) {
        print("Under bound")  
        print(paste0(row_current,";",dim(edgelist)[1]))
    }
    result <- list(graph = edgelist, fitness = NULL, PA = NULL, 
                   type = "directed", fit_time = fit_lifetime, 
                   alpha_time = alpha_lifetime,
                   num_period = num_period,
                   period_length = period_length)
    class(result) <- "PAFit_net"
    return(result)
}