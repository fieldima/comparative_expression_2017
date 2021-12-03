#One big script with everything

#fitContinuous on comparative expression data

#Load packages
library(geiger)
library(arbutus)
library(treeio)
library(ggtree)
library(tidyverse)

#Load treedata
calibrated <- readRDS("export/calibrated_trees") %>% discard(~ Nnode2(.x) < 9)

rm_neg <- function (edge){
  ifelse(edge < 0, 1e-06, edge)
}

#Make a function to iterate through each list
adequacy <- function (tree) {
  result <- vector("list", length = length(tree))
  count <- 1
  for (tr in tree){
    tr@phylo$edge.length <- map_dbl(tr@phylo$edge.length, rm_neg)
    phylo <- tr@phylo
    dat <- tr@data %>% select(label, Mean) %>% #testing to see what happens if we just use mean
      slice(1:((nrow(tr@data)+1)/2)) %>% column_to_rownames(var = "label") %>%
      as.matrix()
    dat[is.na(dat)] = 0.01
    fitBM <- fitContinuous(phylo, dat = dat, model = "BM")
    fitOU <- fitContinuous(phylo, dat = dat, model = "OU")
    fitEB <- fitContinuous(phylo, dat = dat, model = "EB")
    aic <- c(fitBM$opt[["aic"]], fitOU$opt[["aic"]], fitEB$opt[["aic"]])
    fit <- ifelse(min(aic) == aic[1], list(c(fitBM, model = "BM")), 
                  ifelse(min(aic) == aic[2], list(c(fitOU, model = "OU")), 
                         list(c(fitEB, model = "EB"))))
    result[[count]] <- fit
    count <- count + 1
  }
  result
}

toy <- calibrated[1]
adequacy(toy)

first <- calibrated[1:2777]
second <- calibrated[2778:5554]
third <- calibrated[5555:8333]

#Measure fit
first_sec <- adequacy(first)
saveRDS(first_sec, "arbutus/first_fit")
second_sec <- adequacy(second)
saveRDS(second_sec, "arbutus/second_fit")
third_sec <- adequacy(third)
saveRDS(third_sec, "arbutus/third_fit")


#Load data
data <- c(first_sec, second_sec, third_sec)

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

b <- df %>% pivot_longer(c(OU, BM, EB), names_to = "model")

b %>% ggplot(aes(model, value)) + geom_col()
ggsave(filename = "arbutus/AIC_plot.png")

run_arb <- function (fits){
  arby <- vector("list", length = length(fits))
  count = 1
  for(fit in fits){
    for(f in fit){
      class(f) <- "gfit"
      arby[[count]] <- arbutus(f)
      count = count + 1
    }
  }
  arby
}

arb_result <- run_arb(data)