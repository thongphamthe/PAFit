
if (FALSE) {
  #setwd("tests")
  last.warning <- NULL
  for (i in 1:1) {
      source("test_Newman_Kong.R")
      if (length(last.warning) > 0)
          break  
  }
}