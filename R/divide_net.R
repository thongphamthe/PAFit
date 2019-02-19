# function to generate simulated network  with PA and node fitness changes over time
.divide_net <-
  function(graph,
           num_period         = 50 , # number of periods
           period_length      = 20,
           g                  = 100,
           name_start = "Experiment_",
           directory  = "./") {
    time_stamp  <- as.vector(graph[,3])
    unique_time <- sort(unique(time_stamp))
    end_time    <- 1 + 1:num_period * period_length
    start_time  <- c(1,end_time[-length(end_time)])
    in_node     <- graph[time_stamp < unique_time[length(unique_time)],2]
    in_node     <- in_node[in_node != -1]
    out_node    <- graph[time_stamp < unique_time[length(unique_time)],1]
    node_id     <- sort(unique(c(in_node,out_node)))
    for (ii in 1:length(start_time)) {
        data_part  <- graph[time_stamp <= unique_time[end_time[ii]],]
        data_part[data_part[,3] <= unique_time[start_time[ii]],3] <- unique_time[start_time[ii]]
        net_part   <- as.PAFit_net(data_part)
        #stat_temp  <- get_statistics(net_part)
        #g          <- max(round(sqrt(stat_temp$deg_max)),50)
        stats_part <- get_statistics(net_part,g = g)
        part_num   <- ii
        save(data_part,net_part,stats_part,part_num,
             num_period,period_length,node_id, file = paste0(directory,name_start,ii,".Rdata"))
    }
}