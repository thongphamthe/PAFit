
generate_simulated_data_from_estimated_model <- function(net_object, net_stat, result, M = 5) {

oopts <- options(scipen = 999)
on.exit(options(oopts))

if (!is(net_object,"PAFit_net"))
  stop("Error: net_object should be a PAFit_net object.")

if (!is(net_stat,"PAFit_data"))
  stop("Error: net_stat should be a PAFit_data object.")

if (!is(result,"Full_PAFit_result"))
  stop("Error: result should be a Full_PAFit_result object.")

M <- as.integer(M) 

if (M < 3)
  stop("Error: M is too small.")  



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

start_deg         <- 0
deg_new           <- rep(0,length(node_id))
names(deg_new)    <- as.numeric(node_id)
deg_new[as.character(as.numeric(labels(deg)[[1]]))] <- deg
deg               <- deg_new
final_deg         <- deg
deg.max           <- as.numeric(max(deg))
unique_time       <- sort(unique(time_stamp))
T                 <- length(unique_time)
N                 <- length(node_id)


existing_node          <- vector(mode = "list",length = T)
new_node_list          <- vector(mode = "list",length = T) # the new nodes appear at that time-step
const_graph_list       <- vector(mode = "list",length = T) #list of non-simulated edges or new nodes at each time-step
edge_from_new_node     <- vector(mode = "list",length = T)  # list of the fixed source node of new edges needed to draw
num_of_new_edges_fixed <- rep(0,length = T) # total number of new edges needed to draw with fixed source node : only in directed case
num_of_new_edges_free  <- rep(0,length = T) # total number of new edges needed to draw free both end


for (i in 1:T) {
  current_graph <- net[net[,3] == unique_time[i],,drop = FALSE]
  in_node_temp  <- current_graph[, 2, drop = FALSE]
  out_node_temp <- current_graph[, 1, drop = FALSE] 
  
  if (i == 1) {
    const_graph_list[[i]] <- current_graph
  } else {
    if (net_type[1] == "undirected") {
      const_index             <- !(in_node_temp %in% existing_node[[i - 1]]) | !(out_node_temp %in% existing_node[[i - 1]])
      const_graph_list[[i]]   <- current_graph[const_index,,drop = FALSE]
      edge_from_new_node[[i]] <- NULL
      num_of_new_edges_fixed[i] <- 0
      num_of_new_edges_free[i]         <- sum(net[,3] == unique_time[i]) - sum(const_index)
    } else if (net_type[1] == "directed") {
      # destination node are new
      const_index             <- !(in_node_temp %in% existing_node[[i - 1]])
      const_graph_list[[i]]   <- current_graph[const_index,,drop = FALSE]
      # source node index: source node is new but destination node is existing
      
      source_new_index          <- (in_node_temp %in% existing_node[[i - 1]]) & !(out_node_temp %in% existing_node[[i - 1]]) 
      source_node_is_new        <- out_node_temp[source_new_index]
      edge_from_new_node[[i]]   <- source_node_is_new
      num_of_new_edges_fixed[i] <- length(source_node_is_new)
      
      # the remaining new edges that have both source and destination nodes as existing nodes
      num_of_new_edges_free[i]      <- sum(net[,3] == unique_time[i]) - sum(const_index)  - sum(source_new_index)
      
      #check:
      if (num_of_new_edges_free[i] != sum((in_node_temp %in% existing_node[[i - 1]]) & (out_node_temp %in% existing_node[[i - 1]]))) {
        print("Mismatch in num of edge free") 
        print(i)
        print(num_of_new_edges_free[i])
        print(sum((in_node_temp %in% existing_node[[i - 1]]) & (out_node_temp %in% existing_node[[i - 1]])))
      }
      
    }
  }
  existing_node[[i]] <- as.numeric(sort(union(in_node_temp[in_node_temp !=  -1],
                                              out_node_temp[out_node_temp != - 1])))
  if (i > 1) {
    new_node_list[[i]] <- setdiff(existing_node[[i]],existing_node[[i - 1]])
    existing_node[[i]] <- union(existing_node[[i]],existing_node[[i - 1]])
    
  }
}



deg_second_max <- max(result$estimate_result$k)

PA <- (0:deg_second_max)^result$estimate_result$alpha
names(PA) <- 0:deg_second_max
PA["0"]   <- 1
PA[as.character(result$estimate_result$k)] <- result$estimate_result$A

f  <- result$estimate_result$f


is_only_PA <- result$estimate_result$only_PA
is_only_f  <- result$estimate_result$only_f

if (net_type[1] == "directed") {is_directed <- TRUE} else {is_directed <- FALSE};

binning_used    <- net_stat$binning
g_used          <- result$estimate_result$g
deg_thresh_used <- result$estimate_result$deg_threshold
p_used          <- result$cv_data$p
stop_cond_used  <- result$estimate_result$stop_cond

graph_list  <- vector(mode = "list", length = M)
stats_list  <- vector(mode = "list", length = M)
result_list <- vector(mode = "list",length = M)

for (mm in 1:M) {
  graph         <- const_graph_list[[1]]  
  in_node_temp  <- graph[, 2, drop = FALSE]
  out_node_temp <- graph[, 1, drop = FALSE] 
  ok_id         <- which(in_node_temp != -1 & out_node_temp != -1)
  deg_vec       <- rep(-1,length(existing_node[[T]]))
  names(deg_vec) <- as.character(existing_node[[T]])
  
  deg_vec[as.character(existing_node[[1]])] <- 0
  
  if (length(ok_id) > 0) {
      if (TRUE == is_directed) {
          temp_vec       <- table(in_node_temp[ok_id])
      } else {temp_vec <- table(c(in_node_temp[ok_id],out_node_temp[ok_id]))}
      deg_vec[names(temp_vec)] <- deg_vec[names(temp_vec)] + temp_vec
  }
  
  for (i in 2:T) {
    
    if (dim(const_graph_list[[i]])[1] > 0) graph <- rbind(graph,const_graph_list[[i]])
    pa_value <- PA[as.character(deg_vec[as.character(existing_node[[i - 1]])])]
    if (sum(is.na(pa_value)) > 0) {
      if (max(deg_vec) <= deg_second_max ) {
        print("Wrong about deg_second_max ")  
      } else pa_value[is.na(pa_value)] <- PA[length(PA)] # only happen when deg_max > observed deg max
    }
    
    fit_value <- f[as.character(existing_node[[i - 1]])]
    
    if (num_of_new_edges_free[i] > 0) {# sampling new edges whose both ends are free
      in_node_new  <- sample(size = num_of_new_edges_free[i], replace = TRUE, x = existing_node[[i - 1]], 
                             prob = pa_value * fit_value /sum(pa_value * fit_value))  
      if (TRUE == is_directed) {
        out_node_new <- sample(size = num_of_new_edges_free[i],replace = TRUE, x = existing_node[[i - 1]]) # sample out nodes uniformly
      } else {
        out_node_new  <- sample(size = num_of_new_edges_free[i], replace = TRUE, x = existing_node[[i - 1]], 
                                prob = pa_value * fit_value /sum(pa_value * fit_value) ) # based on model
      }
      graph        <- rbind(graph,cbind(out_node_new,in_node_new,unique_time[i]))
    }
    
    if (num_of_new_edges_fixed[i] > 0) { # only possible in directed case: sampling the destination nodes
      if (FALSE == is_directed) {print("Wrong in replicating new edge fixed")}
      in_node_new  <- sample(size = num_of_new_edges_fixed[i], replace = TRUE, x = existing_node[[i - 1]], 
                             prob = pa_value * fit_value/sum(pa_value * fit_value))  
      
      graph <- rbind(graph,cbind(edge_from_new_node[[i]],in_node_new,unique_time[i]))
    }
    
    # update degree vector 
    in_node_temp  <- graph[graph[,3] == unique_time[i], 2, drop = FALSE]
    out_node_temp <- graph[graph[,3] == unique_time[i], 1, drop = FALSE] 
    ok_id         <- which(in_node_temp != -1 & out_node_temp != -1)
    if (length(ok_id) > 0) {
        if (TRUE == is_directed) {
          temp_vec       <- table(in_node_temp[ok_id])
        } else {temp_vec <- table(c(in_node_temp[ok_id],out_node_temp[ok_id]))}
    
        if (length(new_node_list[[i]]) > 0) {
           # check these nodes are really new
           if (sum(deg_vec[as.character(new_node_list[[i]])] != -1) > 0) {
               print("Wrong in updating degree vec")   
           }
           deg_vec[as.character(new_node_list[[i]])] <- 0
        }
        deg_vec[names(temp_vec)] <- deg_vec[names(temp_vec)] + temp_vec
    }
  }
  graph_list[[mm]]  <- graph
  one_net <- as.PAFit_net(graph,type = net_type)
  stats_list[[mm]]  <- get_statistics(one_net,binning = binning_used, g = g_used) 
  if (TRUE == is_only_PA) {
    result_list[[mm]] <- only_A_estimate(net_object = one_net,net_stat = stats_list[[mm]],p = p_used, 
                                         stop_cond = stop_cond_used)
  } else if (TRUE == is_only_f) {
    result_list[[mm]] <- only_F_estimate(net_object = one_net,net_stat = stats_list[[mm]],p = p_used, 
                                         stop_cond = stop_cond_used)
  } else {
      result_list[[mm]] <- joint_estimate(net_object = one_net,net_stat = stats_list[[mm]],p = p_used, 
                                      stop_cond = stop_cond_used)
  }
}
supplement_data <- list(existing_node      = existing_node,
                        new_node_list      = new_node_list,
                        const_graph_list   = const_graph_list, 
                        edge_from_new_node = edge_from_new_node,
                        num_of_new_edges_fixed   = num_of_new_edges_fixed,
                        num_of_new_edges_free    = num_of_new_edges_free)

final_result <- list(result_list = result_list, graph_list = graph_list,
                     stats_list = stats_list, supplement_data = supplement_data)
class(final_result) <- "Simulated_Data_From_Fitted_Model"
return(final_result)
}
