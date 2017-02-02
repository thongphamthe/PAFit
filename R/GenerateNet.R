# function to generate simulated network  
GenerateNet <-
function(N, 
         num_seed       = 2      , 
         multiple_node  = 1      , 
         specific_start = NULL   ,
         m              = 1      ,
         prob_m         = FALSE  ,
         increase       = FALSE  , 
         log            = FALSE  , 
         custom_PA      = NULL   ,
         mode           = 1      , 
         alpha          = 1      , 
         beta           = 2      , 
         sat_at         = 100    ,
         offset         = 1      ,
         mode_f         = "gamma", 
         rate           = 0      , 
         shape          = 0      , 
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
   if ((alpha < 0 && mode != 1) || (beta < 0) || (sat_at < 0) || (rate < 0) || (shape < 0) || (m <= 0))
       stop("The parameters must be non-negative")  
   if ((mode[1] != 1) && (mode[1] != 3) && (mode[1] != 2))
       stop("Mode must be 1, 2 or 3 ")  
    graph <- vector("list", N)
    # Node weights in the BA model. 
    switch(mode_f[1], 
    gamma = {
    # gamma distribution
      if (shape*rate > 0) 
          fitness <- rgamma(N,rate = rate,shape = shape)
      else fitness <- rep(1,N)},
     log_normal ={
     # log_normal distribution
        fitness <- rlnorm(N, meanlog = meanlog, sdlog = sdlog)
     },
    power_law = {
        fitness <- rpareto(N, scale = scale_pareto, shape = shape_pareto)
    },{
    stop('mode_f must be either "gamma", "log_normal" or "power_law"')
    }
    )
    names(fitness) <- 1:N
    for (n in 2:num_seed)
         graph[[n]] <- n - 1
    
    degree        <- rep(0,num_seed)
    names(degree) <- 1:num_seed 
    for (i in 1:num_seed) {
        count_degree <-table(graph[[i]][graph[[i]] <= i])
        degree[labels(count_degree)[[1]]] <- degree[labels(count_degree)[[1]]] + 
                                             count_degree
    }
    P         <- degree
    P[P == 0] <- offset
    seed.graph.size <- length(P)
    count <- 0
    # n is the size of the current network
    n <- seed.graph.size
    while (n < N) {	
        # at each time_step: save n
        n_old <- n
        if (!is.null(custom_PA)) {
            A                 <- custom_PA
            final_A           <- A[length(A)]  
            temp              <- A[P + 1]
           #print(temp)
            temp[is.na(temp)] <- final_A
            P.sum             <- sum(temp*fitness[1:n_old])
            node.weights      <- temp*fitness[1:n_old]/P.sum  
        }
        else
        if (mode[1] == 1) {
            P.sum        <- sum(P^alpha*fitness[1:n_old])
            node.weights <- P^alpha*fitness[1:n_old]/P.sum
        } else if (mode[1] == 2) {
            temp         <- pmin(P,sat_at)^alpha
            P.sum        <- sum(temp*fitness[1:n_old])
            node.weights <- temp*fitness[1:n_old]/P.sum
        } else {
            temp         <- alpha*(log(P))^beta + 1 
            P.sum        <- sum(temp*fitness[1:n_old])
            node.weights <- temp*fitness[1:n_old]/P.sum
        }
        for (i in 1:multiple_node){
        #Allow duplication:  
            if (TRUE == increase) {
                count <- count + 1  
                nodes <- sort(sample(1:n_old,size = ifelse(log,max(round(log(count)),1),count),prob = node.weights, replace = TRUE))      
            } else {
            if (prob_m == TRUE) {
                num_edge_temp <- rpois(1,lambda = m)
            if (num_edge_temp > 0)
                nodes <- sort(sample(1:n_old,num_edge_temp,prob = node.weights, replace = TRUE))
            else nodes <- NULL
            }
           else
               nodes <- sort(sample(1:n_old,size = m,prob = node.weights, replace = TRUE)) 
            }
            ##########################################
            degree  <- c(degree,0)
            if (0 != length(nodes)) {
                temp    <- table(nodes)
                graph[[n+1]]  <- c(graph[[n+1]], nodes)
                for(i in 1:length(temp)) { 
                    num_edge            <- as.numeric(temp[i]) 
                    node_name           <- as.numeric(labels(temp[i]))
                    degree[node_name]   <- degree[node_name] + num_edge # Update degrees.
                }
            }
            n <- n + 1
            if (n == N) break  
        }
        P <- degree
        P[degree == 0] <- offset     
    }
    
    
    num_of_edge          <- sum(unlist(lapply(graph,function(x) length(as.vector(x)))))
    edge_list            <- matrix(nrow = num_of_edge,ncol = 3,0)
    sum_m                <- 0
    current_time_step    <- 0
    i                    <- 1
    flag_specific_start  <- !is.null(specific_start)
    break_flag           <- FALSE
    while (i <= N) {
        if (break_flag == TRUE)
            break  
        #print(i)
        if (i <= num_seed) {  
            m_t <- length(graph[[i]][graph[[i]] <= i])
            if (m_t > 0) { 
                temp  <- as.vector(graph[[i]])
                edge_list[(sum_m + 1):(sum_m + m_t),3]   <-  0
                edge_list[(sum_m + 1):(sum_m + m_t),1]   <-  i
                edge_list[(sum_m + 1):(sum_m  +  m_t),2] <- temp 
                sum_m                                    <- sum_m + m_t
            }
           i <- i + 1 
           if (i > N) { break_flag = TRUE;break}  
        }
        else {
            current_time_step <- current_time_step + 1;  
            if (1 == flag_specific_start) {
                for (j in 1:specific_start) {
                    if (i > N) {break_flag = TRUE;break}  
                    m_t <- length(graph[[i]][graph[[i]] <= i])
                    if (m_t > 0) { 
                      temp  <- as.vector(graph[[i]])
                      edge_list[(sum_m + 1):(sum_m + m_t),3]   <-  current_time_step - 1
                      edge_list[(sum_m + 1):(sum_m + m_t),1]   <-  i
                      edge_list[(sum_m + 1):(sum_m  +  m_t),2] <- temp 
                      sum_m                                    <- sum_m + m_t
                    }
                    i <- i + 1  
                    if (i > N) { break_flag = TRUE;break}  
                }
              flag_specfific_start <- 0;
            } else 
            for (j in 1:multiple_node) {
              if (i > N) {break_flag = TRUE; break}    
              m_t <- length(graph[[i]][graph[[i]] <= i])
                  if (m_t > 0) { 
                      temp  <- as.vector(graph[[i]])
                      edge_list[(sum_m + 1):(sum_m + m_t),3]   <-  current_time_step
                      edge_list[(sum_m + 1):(sum_m + m_t),1]   <-  i
                      edge_list[(sum_m + 1):(sum_m  +  m_t),2] <-  temp 
                      sum_m                                    <- sum_m + m_t
                  }
            i <- i + 1
            if (i > N) { break_flag = TRUE;break}  
           }  
      }
  }
  return(list(graph = edge_list, fitness = fitness))
}
