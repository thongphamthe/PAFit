# function to summarize statistics from a growing network 
get_statistics <-
function(net_object ,
         only_PA       = FALSE , only_true_deg_matrix = FALSE,
         binning       = TRUE  , g                    = 50   ,
         deg_threshold = 0     , 
         compress_mode = 0     , compress_ratio       = 0.5  , 
         custom_time   = NULL){
    #net               <- as.matrix(net)
  oopts <- options(scipen = 999)
  on.exit(options(oopts))
  
    if (!is(net_object,"PAFit_net"))
        stop("Error: net_object should be a PAFit_net object.")
    net               <- net_object$graph
    net_type          <- net_object$type
    net               <- net[order(net[,3], decreasing = FALSE),]
    time_stamp        <- as.vector(net[,3])
    in_node           <- as.vector(net[,2])
    out_node          <- as.vector(net[,1])
    out_node          <- out_node
    node_id           <- as.numeric(sort(union(in_node[in_node !=  -1],out_node[out_node != - 1])))
    
    ok_id <- which(in_node != -1 & out_node != -1)
    if (net_type[1] == "directed") {
        deg           <- table(in_node[ok_id])
    } else
    if (net_type[1] == "undirected")
        deg           <- table(c(in_node[ok_id],out_node[ok_id]))     
        #deg           <- table(as.vector(as.matrix(net[,1:2])))        
    start_deg         <- 0
    deg_new           <- rep(0,length(node_id))
    names(deg_new)    <- as.numeric(node_id)
    deg_new[as.character(as.numeric(labels(deg)[[1]]))] <- deg
    deg               <- deg_new
    final_deg         <- deg
    deg.max           <- as.numeric(max(deg))
    if (start_deg >= deg.max)
        stop("Starting degree too large!") else 
    if (start_deg < 0)
        stop("Negative starting degree!")
    unique_time       <- sort(unique(time_stamp))
    T                 <- length(unique_time)
    N                 <- length(node_id)
    if (only_true_deg_matrix == TRUE)
        binning <- FALSE  
    ##############  binning #########################
    #We have to cover from start_deg degree to deg.max degree, that is an interval with length deg.max - start_deg + 1
    if ((TRUE == binning) && (g > 0) && (g <= deg.max - start_deg + 1)) {
        if (1 == g) {
            base            <- deg.max - start_deg + 1
            interval_length <- deg.max - start_deg + 1
        } else {
            #find the base that gives exactly g bin
            is.warn <- options()$warn 
            options(warn = -1) #temporily supress warning
            ff <- function(x){deg.max - start_deg + 1.0 - sum(floor(x^(0:(g - 1))))}
            base      <- uniroot(ff,interval = c(1 + 1e-15,deg.max - start_deg + g + 1.1),tol = .Machine$double.eps)$root
            options(warn = is.warn)
            interval_length <- floor(base^(0:(g-1)))
        }
    } else if ((FALSE == binning) || (0 == g) || (g > deg.max - start_deg + 1)) {
        g               <- deg.max - start_deg + 1
        interval_length <- rep(1,g)
        base            <- 1
    }
    
    bin_vector   <-rep(g + start_deg - 1, deg.max + 1)  # degree 0 to deg.max
    # The right-end degree of the bins
    begin_deg   <- c(start_deg,start_deg + cumsum(interval_length)[-g])
    end_deg     <- begin_deg + interval_length - 1
    if (start_deg > 0)
        bin_vector[1:start_deg]  <- 0:(start_deg - 1)
    for (i in 1:g) 
        bin_vector[(begin_deg[i]:end_deg[i]) + 1]  <- i + start_deg - 1
    ########### Compress Time stamp #########################
    
    if (1 == compress_mode[1]) {
        T_compressed           <- round(compress_ratio*(T - 1))
        compressed_unique_time <- unique_time[floor(seq(1,T - 1,length.out = T_compressed))]
    } else if (2 == compress_mode[1]){
        edge_cumsum             <- cumsum(as.vector(table(time_stamp))) 
        edge_ratio              <- edge_cumsum/edge_cumsum[T]
        compressed_unique_time  <- unique_time[which(edge_ratio >= 1 - compress_ratio)]
        temp                    <- edge_ratio[which(edge_ratio < 1 - compress_ratio)]
        compress_ratio           <- temp[length(temp)]
        T_compressed            <- length(compressed_unique_time)
    }  else if (3 == compress_mode[1]) {
        compressed_unique_time <- sort(unique(custom_time))
        T_compressed           <- length(compressed_unique_time)
        compress_ratio          <- length(compressed_unique_time)/(T - 1)
    } else {
        #No Time compression
        compress_ratio   <- 1
        T_compressed    <- T
        compressed_unique_time <- unique_time
    }
    if (net_type[1] == "directed")
        first_edge       <- table(in_node[ok_id][time_stamp[ok_id] == unique_time[1]]) else 
    if (net_type[1] == "undirected") {
        net_temp   <- as.matrix(net[ok_id,][time_stamp[ok_id] == unique_time[1],1:2])
        #net_temp   <- net_temp[ok_id,1:2]
        first_edge <- table(as.vector(net_temp))
    }
       
    first_deg        <- rep(0,N)
    increase         <- rep(0,N)
    names(increase)  <- as.numeric(node_id)
    names(first_deg) <- as.numeric(node_id)
    first_deg[as.character(as.numeric(labels(first_edge)[[1]]))] <- first_edge
    # this is not the true number of new edges, due to the first appearance of a node together with some edges
    # but this can be used as an crude first step selection
    # the accurate selection can be done after
    inc              <- deg - first_deg  
    increase         <- inc
    initial_nodes    <- length(first_edge)
    if (FALSE == only_PA) {
        pos_temp          <- inc >= deg_threshold
        if (sum(pos_temp) == 0)
            stop("Degree threshold is too high. Please decrease degree threshold.")  
        f_position        <- node_id[pos_temp]
    } else f_position        <- NULL
    
    node_id_old <- node_id
    if (FALSE == only_PA)
        node_id     <- node_id[inc >= deg_threshold]
    
    N_new            <- length(node_id) 
    degree_appear    <- rep(0,deg.max + 1)
    if (only_true_deg_matrix == TRUE) {
        Sum_m_k          <- vector() 
        n_tk             <- matrix(0,0,0)
        m_tk             <- matrix(0,0,0)
        m_t              <- vector() 
    } else {
        Sum_m_k          <- rep(0,start_deg + g)
        n_tk             <- matrix(0,nrow = T_compressed - 1, ncol = start_deg + g)
        m_tk             <- matrix(0,nrow = T_compressed - 1, ncol = start_deg + g)
        m_t              <- rep(0,T_compressed - 1)  
    }

    if (FALSE == only_PA) {
        if (only_true_deg_matrix == FALSE)  {
            offset_tk               <- matrix(0,nrow = T_compressed - 1, ncol = start_deg + g) 
            offset_m_tk             <- matrix(0,nrow = T_compressed - 1, ncol = start_deg + g) 
            z_j                     <- rep(0,N_new)
            appear_time             <- rep(0,N_new) #index of m_t, start from 1
        } else {
          offset_tk               <- matrix(0,0,0)
          offset_m_tk             <- matrix(0,0,0) 
          z_j                     <- vector()    
          appear_time             <- vector() #index of m_t, start from 1
        }
        node_degree <- matrix(-1,nrow = T_compressed - 1, ncol = N_new) 
    } else {
        offset_tk               <- matrix(0,0,0)
        offset_m_tk             <- matrix(0,0,0) 
        z_j                     <- vector()  
        node_degree             <- matrix(0,0,0)
        appear_time             <- vector() #index of m_t, start from 1
    }
    if (net_type[1] == "directed")
        undirected  = 0
    else undirected = 1;
    max_node_id = max(node_id_old);
    only_PA_num = ifelse(only_PA,1,0);
    only_true_deg_matrix_num = ifelse(only_true_deg_matrix,1,0)
    #print(time_stamp)
    #print(unique_time)
    center_k <- rep(0,g)
    
    #print(bin_vector)
    
    #print(length(bin_vector))
    
    #print(dim(node_degree))
    #print(length(Sum_m_k))
    #print(dim(m_tk))
    #print(dim(n_tk))
    #print(max_node_id)
    #print(length(node_id_old))
    
    
    #print(system.time({
    #for (j in 1:length(new_id)) {
    #    in_node[in_node == node_id_old[j]]   <- j
    #    out_node[out_node == node_id_old[j]] <- j
    #}
    #  .my_replace(in_node,out_node,node_id_old)
    #}))
    
    #node_id_old <- new_id

    
    .get_stats(time_stamp,unique_time, in_node, out_node, node_id_old, node_id,bin_vector, max_node_id, undirected, 
              only_PA_num,              
              compressed_unique_time,
              Sum_m_k,n_tk,m_tk,m_t,offset_tk,z_j,node_degree,offset_m_tk,only_true_deg_matrix_num, deg.max, center_k,
              appear_time)
    
    center_k[center_k == 0] <- begin_deg[center_k == 0]
    #if (center_k[length(center_k)] == 0)
    #    center_k[length(center_k)] <- center_k[length(center_k) - 1]  
    
    if (FALSE == only_PA) {
      if (only_true_deg_matrix == FALSE) {   
          names(z_j)            <- as.numeric(node_id)
          #print(only_true_deg_matrix)
          #print("start:")
          #print(z_j)
          names(appear_time)    <- as.numeric(node_id)
          #print(appear_time)
          #print("Stop!")
      }
      colnames(node_degree) <- as.numeric(node_id)
      
    }
   
    names(node_id)    <- as.numeric(node_id)
    #now perform the final selection
    true                           <- which(z_j >= deg_threshold)
    
    if (FALSE == only_PA) {
        if (only_true_deg_matrix == FALSE) {  
            if (length(true) == 0)
                stop("Degree threshold is too high. Please decrease degree threshold.")  
            increase[inc >= deg_threshold] <- z_j
            z_j                            <- z_j[true]
            f_position                     <- f_position[true]
        }
    }
    if (only_true_deg_matrix == FALSE) 
        node_degree                    <- node_degree[,true,drop = FALSE]
    

    result  <- list(offset_tk = offset_tk, offset_m_tk = offset_m_tk, net_type = net_type[1], 
                    n_tk = n_tk,m_tk = m_tk, bin_vector = bin_vector, center_k = center_k, 
                    sum_m_k = Sum_m_k,
                    node_degree = node_degree,m_t = m_t,z_j = z_j, initial_nodes = initial_nodes,
                    deg_thresh = deg_threshold, final_deg = final_deg, only_PA = only_PA, 
                    increase = increase, start_deg = start_deg, 
                    binning = binning, g = g, 
                    compress_mode = compress_mode[1], 
                    f_position = as.numeric(f_position), 
                    compressed_unique_time = compressed_unique_time, 
                    begin_deg = begin_deg, end_deg = end_deg,
                    interval_length = interval_length,
                    node_id = as.numeric(node_id_old), N = N, T = T, 
                    t_compressed = T_compressed,
                    deg_max = deg.max, compress_ratio = compress_ratio , 
                    custom_time = custom_time, 
                    only_true_deg_matrix = only_true_deg_matrix,
                    appear_time = appear_time)
    class(result) <- "PAFit_data"
   
    return(result)
}

.onUnload <- function (libpath) {
  library.dynam.unload("PAFit", libpath)
}
