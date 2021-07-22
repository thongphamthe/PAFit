library(PAFit)

testNet1 <- generate_BA(N = 1000, m = 1, alpha = 0.75)
result   <- PAFit_oneshot(testNet1)
result$final_result$regress_res
