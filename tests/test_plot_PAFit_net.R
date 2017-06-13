

library("PAFit")
net <- generate_BB(N = 100, multiple_node = 20, m = 1)

plot(net)

plot(net, slice = 3)


u <- as.PAFit_net(coauthor.net, type = "undirected")
plot(u)
plot(u, slice = 10)
#plot(network.extract(b,at = 6))

plot(net, plot = "PA")
plot(net, plot = "fit")
