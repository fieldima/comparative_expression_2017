#Make histogram for AICs

#libraries
library(tidyverse)

#Load data
toy_test <- readRDS("arbutus/toy_fit")

#First need to pull out model names
test <- unlist(toy_test[[1]]) 

model_count <- function (fit) {
  ou = 0
  bm = 0
  eb = 0
  for(f in fit){
    vec <- unlist(f)
    ifelse(vec$model == "OU", ou <- ou + 1, ifelse(vec$model == "BM", bm <- bm + 1, eb <- eb + 1))
  }
  data.frame(OU = ou, BM = bm, EB = eb)
}

a <- model_count(toy_test)

b <- a %>% pivot_longer(c(OU, BM, EB), names_to = "model")

b %>% ggplot(aes(model, value)) + geom_col()
