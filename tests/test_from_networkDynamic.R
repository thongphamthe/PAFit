library(PAFit)


net <- generate_BA(N = 100, multiple_node = 20, m = 1)
summary(net)
system.time(b <- to_networkDynamic(net))

back <- from_networkDynamic(b)

