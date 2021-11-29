#run Arbutus

library(arbutus)
library(geiger)

#Load test data
toy_test <- readRDS("arbutus/toy_fit")

run_arb <- function (fits){
arb <- vector("list", length = length(toy_test))
count = 1
for(toy in toy_test){
  for(t in toy){
    class(t) <- "gfit"
    arb[[count]] <- arbutus(t)
    count = count + 1
  }
}
arb
}

run1 <- run_arb(toy_test)
