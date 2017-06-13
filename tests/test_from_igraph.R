
library(PAFit)

net <- generate_BA(N = 100, multiple_node = 20, m = 1)

system.time(b <- to_igraph(net))


back <- from_igraph(b)

#plot(network.extract(b,at = 6))
