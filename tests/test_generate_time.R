if (FALSE) {
  set.seed(1)
  library("PAFit")
   net        <- .generate_net_time(num_seed          = 10    , # number of nodes at the initial stage
                                   num_period         = 50    , # number of periods
                                   period_length      = 10    , # number of time-steps inside each period
                                   multiple_node      = 1     , # number of new nodes appears at each time-steps
                                   m                  = 1     , # number of new edges per new nodes at each time-step
                                   alpha_start        = 0.25  , # starting value of alpha at the start of the period
                                   alpha_end          = 0.75  , # ending value of alpha at the end of the period
                                   s                  = 10    , # inverse variance of gamma where new node fitnesses will be drawn from
                                   h                  = 0.15)
  
}
