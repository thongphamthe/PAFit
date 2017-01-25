#Fit Multinomial
FitMultinomial <- function(true,dat){
  true[true == 0] <- 1
  return(sum(dat*log(true)))
}