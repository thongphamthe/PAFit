library(PAFit)

mean_array <- rep(0,10)
for (i in 1:10) {
testNet1 <- generate_BA(N = 500, m = 1, alpha = 0.75)
result   <- PAFit_oneshot(testNet1)
mean_array[i] <- result$alpha
}

mean(mean_array)
sd(mean_array)/sqrt(10)
