#Make histogram for AICs

#libraries
library(tidyverse)

#Load data
data <- c(readRDS("arbutus/first_fit"),readRDS("arbutus/second_fit"),readRDS("arbutus/third_fit"))

#First need to pull out model names
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

df <- model_count(data)
saveRDS(data, file = "arbutus/fitdata")

b <- df %>% pivot_longer(c(df.OU.2, df.BM.2, df.EB.2), names_to = "model")

b %>% ggplot(aes(model, value)) + geom_col()
