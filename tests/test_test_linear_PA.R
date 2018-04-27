
library("PAFit")
set.seed(1)

for (iii in 1:5) {
    net   <- generate_BA(N = 100)
    stats <- get_statistics(net, only_PA = TRUE)
    u     <- test_linear_PA(stats$final_deg)
    print(u)
}
