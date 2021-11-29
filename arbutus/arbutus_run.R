#run Arbutus

library(arbutus)
library(geiger)

#Load test data
toy_test <- readRDS("arbutus/toy_fit")

for(toy in toy_test){
  for(t in toy){
    arbutus(t)
  }
}

toy <- toy_test[[1]]
t <- toy[[1]]
arbutus(t[-5])
gfit(t)
